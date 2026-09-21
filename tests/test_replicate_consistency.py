"""凝胶重复实验一致性评估：全局条带聚类、共识谱、逐泳道标记与稳定性评分的
单元 + API 端到端测试。"""
from fastapi.testclient import TestClient

from app.main import app
from app.replicate import assess_replicate_consistency
from app.models import ReplicateConsistencyRequest

client = TestClient(app)


# ---------------------------------------------------------------- 构造工具

def _req(lanes, tolerance=10, min_occurrences=2):
    return ReplicateConsistencyRequest(
        lanes=[{"name": n, "observed_bands": b} for n, b in lanes],
        band_match_tolerance=tolerance,
        min_occurrences=min_occurrences,
    )


def _lane_map(res):
    return {item.lane: item for item in res.lanes}


# ---------------------------------------------------------------- 纯算法

def test_replicate_basic_consensus_and_scores():
    # 三泳道重复：C1 {98,100,102} 中位 100、C2 {200,200,210} 中位 200、
    # C3 {300,310} 中位 305，出现次数均 ≥ 2 全部入共识谱
    res = assess_replicate_consistency(_req([
        ("lane-1", [100, 200, 300]),
        ("lane-2", [102, 200, 310]),
        ("lane-3", [98, 210]),
    ]))
    assert res.lane_count == 3
    assert res.total_band_count == 8
    assert res.cluster_count == 3
    assert res.consensus_band_count == 3
    assert res.low_frequency_cluster_count == 0
    assert res.consensus_profile == [100.0, 200.0, 305.0]

    c1, c2, c3 = res.clusters
    assert c1.cluster_id == 1 and c1.median_length == 100.0
    assert c1.member_count == 3 and c1.occurrence_count == 3
    assert c1.in_consensus and c1.lanes == ["lane-1", "lane-2", "lane-3"]
    # 簇内成员按（长度、泳道名）升序
    assert [(m.length, m.lane) for m in c1.members] == [
        (98.0, "lane-3"), (100.0, "lane-1"), (102.0, "lane-2"),
    ]
    assert (c1.min_length, c1.max_length) == (98.0, 102.0)
    assert c3.median_length == 305.0 and c3.occurrence_count == 2

    lanes = _lane_map(res)
    # lane-1：全部命中，偏差 |100-100|+|200-200|+|300-305| = 5
    assert lanes["lane-1"].matched_consensus_count == 3
    assert lanes["lane-1"].anomaly_count == 0
    assert lanes["lane-1"].consensus_deviation == 5.0
    assert lanes["lane-1"].mean_deviation == 1.67
    assert lanes["lane-1"].stability_score == 37.45
    # lane-2：偏差 2+0+5 = 7
    assert lanes["lane-2"].consensus_deviation == 7.0
    assert lanes["lane-2"].mean_deviation == 2.33
    assert lanes["lane-2"].stability_score == 30.03
    # lane-3：缺失共识条带 305，偏差 2+10 = 12
    assert lanes["lane-3"].missing_bands == [305.0]
    assert lanes["lane-3"].missing_count == 1
    assert lanes["lane-3"].anomaly_count == 1
    assert lanes["lane-3"].stability_score == 9.52
    # 整体评分为各泳道评分（响应值）的均值
    assert res.mean_stability_score == round((37.45 + 30.03 + 9.52) / 3, 2)
    assert res.score_reason
    for item in res.lanes:
        assert item.description
        assert item.anomaly_count == (
            item.missing_count + item.extra_count + item.low_frequency_count
        )


def test_replicate_extra_bands_are_surplus_within_consensus_cluster():
    # lane-1 在共识簇 {100,100,105} 中贡献两条：代表 100（距中位最近），105 为额外
    res = assess_replicate_consistency(_req([
        ("lane-1", [100, 105, 200]),
        ("lane-2", [100, 200]),
    ]))
    assert res.consensus_profile == [100.0, 200.0]
    lanes = _lane_map(res)
    assert lanes["lane-1"].extra_bands == [105.0]
    assert lanes["lane-1"].extra_count == 1
    assert lanes["lane-1"].consensus_deviation == 0.0
    # 评分 = 100 × 2/(2+1) ÷ (1+0)
    assert lanes["lane-1"].stability_score == round(100 * 2 / 3, 2)
    assert lanes["lane-2"].extra_bands == []
    assert lanes["lane-2"].stability_score == 100.0


def test_replicate_low_frequency_clusters_and_bands():
    # 900、950 各自只出现在一条泳道（出现次数 1 < 2）-> 低频簇
    res = assess_replicate_consistency(_req([
        ("lane-1", [100, 900]),
        ("lane-2", [100]),
        ("lane-3", [100, 950]),
    ], tolerance=5))
    assert res.cluster_count == 3
    assert res.consensus_band_count == 1
    assert res.low_frequency_cluster_count == 2
    assert res.consensus_profile == [100.0]
    low_clusters = [c for c in res.clusters if not c.in_consensus]
    assert [c.median_length for c in low_clusters] == [900.0, 950.0]
    assert all(c.occurrence_count == 1 for c in low_clusters)

    lanes = _lane_map(res)
    assert lanes["lane-1"].low_frequency_bands == [900.0]
    assert lanes["lane-3"].low_frequency_bands == [950.0]
    assert lanes["lane-2"].low_frequency_bands == []
    # 评分 = 100 × 1/(1+0+1) ÷ (1+0) = 50
    assert lanes["lane-1"].stability_score == 50.0
    assert lanes["lane-2"].stability_score == 100.0


def test_replicate_median_of_even_sized_cluster():
    # 偶数条簇取中间两条均值：(100+110)/2 = 105
    res = assess_replicate_consistency(_req([
        ("lane-1", [100]),
        ("lane-2", [110]),
    ], tolerance=10, min_occurrences=2))
    assert res.consensus_profile == [105.0]
    lanes = _lane_map(res)
    assert lanes["lane-1"].consensus_deviation == 5.0
    assert lanes["lane-2"].consensus_deviation == 5.0


def test_replicate_single_linkage_chains_into_one_cluster():
    # 链式可传递：100~105~110 相邻差 ≤ 5 -> 同一簇，中位 105
    res = assess_replicate_consistency(_req([
        ("lane-1", [100]),
        ("lane-2", [105]),
        ("lane-3", [110]),
    ], tolerance=5, min_occurrences=2))
    assert res.cluster_count == 1
    assert res.clusters[0].median_length == 105.0
    assert res.clusters[0].occurrence_count == 3
    assert res.consensus_profile == [105.0]


def test_replicate_min_occurrences_one_makes_every_cluster_consensus():
    res = assess_replicate_consistency(_req([
        ("lane-1", [100, 900]),
        ("lane-2", [100]),
    ], tolerance=5, min_occurrences=1))
    assert res.low_frequency_cluster_count == 0
    assert res.consensus_profile == [100.0, 900.0]
    lanes = _lane_map(res)
    # 900 的簇只有 lane-1 贡献 -> lane-2 缺失
    assert lanes["lane-2"].missing_bands == [900.0]
    assert lanes["lane-1"].missing_bands == []


def test_replicate_no_consensus_when_all_clusters_below_threshold():
    # 两条泳道条带完全不同且 min_occurrences=2：全部簇低频，共识谱为空
    res = assess_replicate_consistency(_req([
        ("lane-1", [100]),
        ("lane-2", [200]),
    ], tolerance=5, min_occurrences=2))
    assert res.consensus_profile == []
    assert res.consensus_band_count == 0
    assert res.low_frequency_cluster_count == 2
    lanes = _lane_map(res)
    assert lanes["lane-1"].low_frequency_bands == [100.0]
    assert lanes["lane-1"].matched_consensus_count == 0
    assert lanes["lane-1"].mean_deviation == 0.0
    # 无命中 -> 评分 0
    assert lanes["lane-1"].stability_score == 0.0
    assert res.mean_stability_score == 0.0


def test_replicate_tie_break_by_length_then_lane_rank_is_canonical():
    # 等长条带跨泳道并列：成员按泳道名（规范名次）排列，与提交顺序无关
    res = assess_replicate_consistency(_req([
        ("zeta", [100]),
        ("alpha", [100]),
        ("mid", [100]),
    ], tolerance=0, min_occurrences=2))
    assert res.cluster_count == 1
    assert [m.lane for m in res.clusters[0].members] == ["alpha", "mid", "zeta"]
    assert [item.lane for item in res.lanes] == ["alpha", "mid", "zeta"]
    assert all(item.stability_score == 100.0 for item in res.lanes)


def test_replicate_result_is_input_order_independent():
    req_a = _req([
        ("lane-1", [100, 200, 300]),
        ("lane-2", [102, 200, 310]),
        ("lane-3", [98, 210]),
    ])
    req_b = _req([
        ("lane-3", [210, 98]),
        ("lane-2", [310, 102, 200]),
        ("lane-1", [300, 100, 200]),
    ])
    assert assess_replicate_consistency(req_a).model_dump() == \
        assess_replicate_consistency(req_b).model_dump()


def test_replicate_perfect_lanes_score_100():
    res = assess_replicate_consistency(_req([
        ("lane-1", [100, 200]),
        ("lane-2", [100, 200]),
        ("lane-3", [100, 200]),
    ], tolerance=5, min_occurrences=2))
    assert all(item.stability_score == 100.0 for item in res.lanes)
    assert res.mean_stability_score == 100.0
    assert all(item.anomaly_count == 0 for item in res.lanes)


def test_replicate_counts_are_consistent():
    res = assess_replicate_consistency(_req([
        ("lane-1", [100, 105, 900]),
        ("lane-2", [100, 200]),
        ("lane-3", [98, 200, 950]),
    ], tolerance=5, min_occurrences=2))
    for item in res.lanes:
        # 每条条带恰好归入一类：代表命中 / 额外 / 低频
        assert (
            item.matched_consensus_count + item.extra_count
            + item.low_frequency_count == item.band_count
        )
        assert item.matched_consensus_count + item.missing_count == \
            res.consensus_band_count
        assert item.missing_count == len(item.missing_bands)
        assert item.extra_count == len(item.extra_bands)
        assert item.low_frequency_count == len(item.low_frequency_bands)
    # 各簇出现次数 = 贡献泳道数，成员数 = 各泳道条带数之和
    for cluster in res.clusters:
        assert cluster.occurrence_count == len(cluster.lanes)
        assert cluster.member_count == len(cluster.members)
    assert sum(c.member_count for c in res.clusters) == res.total_band_count


# ---------------------------------------------------------------- API 端到端

def _payload(**kw):
    payload = {
        "lanes": [
            {"name": "lane-1", "observed_bands": [100, 200, 300]},
            {"name": "lane-2", "observed_bands": [102, 200, 310]},
            {"name": "lane-3", "observed_bands": [98, 210]},
        ],
        "band_match_tolerance": 10,
        "min_occurrences": 2,
    }
    payload.update(kw)
    return payload


def test_api_replicate_consistency_full_response():
    res = client.post("/api/v1/digest/replicate-consistency", json=_payload())
    assert res.status_code == 200
    body = res.json()

    assert body["band_match_tolerance"] == 10
    assert body["min_occurrences"] == 2
    assert body["lane_count"] == 3
    assert body["total_band_count"] == 8
    assert body["cluster_count"] == 3
    assert body["consensus_band_count"] == 3
    assert body["low_frequency_cluster_count"] == 0
    assert body["consensus_profile"] == [100.0, 200.0, 305.0]

    assert [c["cluster_id"] for c in body["clusters"]] == [1, 2, 3]
    first = body["clusters"][0]
    assert first["median_length"] == 100.0
    assert first["occurrence_count"] == 3
    assert first["in_consensus"] is True
    assert first["lanes"] == ["lane-1", "lane-2", "lane-3"]
    assert [(m["length"], m["lane"]) for m in first["members"]] == [
        (98.0, "lane-3"), (100.0, "lane-1"), (102.0, "lane-2"),
    ]

    assert [l["lane"] for l in body["lanes"]] == ["lane-1", "lane-2", "lane-3"]
    lane3 = body["lanes"][2]
    assert lane3["missing_bands"] == [305.0]
    assert lane3["anomaly_count"] == 1
    assert lane3["stability_score"] == 9.52
    assert body["mean_stability_score"] == 25.67
    assert body["score_reason"]


def test_api_replicate_consistency_order_independent():
    body_a = client.post(
        "/api/v1/digest/replicate-consistency", json=_payload(),
    ).json()
    payload = _payload(lanes=[
        {"name": "lane-3", "observed_bands": [210, 98]},
        {"name": "lane-1", "observed_bands": [300, 100, 200]},
        {"name": "lane-2", "observed_bands": [102, 200, 310]},
    ])
    body_b = client.post(
        "/api/v1/digest/replicate-consistency", json=payload,
    ).json()
    assert body_a == body_b


# ---------------------------------------------------------------- 422 校验

def _expect_422(payload, code_part=None):
    res = client.post("/api/v1/digest/replicate-consistency", json=payload)
    assert res.status_code == 422
    body = res.json()
    assert body["code"] == "validation_failed"
    assert body["errors"]
    if code_part:
        assert any(code_part in e["type"] for e in body["errors"])
    return body


def test_replicate_422_duplicate_and_blank_lane_names():
    _expect_422(_payload(lanes=[
        {"name": "lane-1", "observed_bands": [100]},
        {"name": "lane-1", "observed_bands": [200]},
    ]), "duplicate_lane_name")
    _expect_422(_payload(lanes=[
        {"name": "  ", "observed_bands": [100]},
        {"name": "lane-2", "observed_bands": [200]},
    ]), "invalid_lane_name")


def test_replicate_422_empty_lists_and_missing_fields():
    _expect_422(_payload(lanes=[]))
    _expect_422(_payload(lanes=[{"name": "lane-1", "observed_bands": []}]))
    payload = _payload()
    del payload["min_occurrences"]
    _expect_422(payload)
    del payload["band_match_tolerance"]
    _expect_422(payload)
    _expect_422({})


def test_replicate_422_negative_values_and_bad_types():
    _expect_422(_payload(band_match_tolerance=-1))
    _expect_422(_payload(min_occurrences=0))
    _expect_422(_payload(min_occurrences=-2))
    _expect_422(_payload(min_occurrences="two"))
    _expect_422(_payload(lanes=[
        {"name": "lane-1", "observed_bands": [100, -5]},
    ]))
    _expect_422(_payload(band_match_tolerance="wide"))


def test_replicate_422_min_occurrences_exceeds_lane_count():
    _expect_422(_payload(min_occurrences=4), "min_occurrences_too_large")
    # 等于泳道数时合法
    res = client.post(
        "/api/v1/digest/replicate-consistency", json=_payload(min_occurrences=3),
    )
    assert res.status_code == 200


def test_replicate_422_too_many_lanes_or_bands(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_REPLICATE_LANES", 3)
    _expect_422(_payload(lanes=[
        {"name": f"lane-{i}", "observed_bands": [100 + i]} for i in range(4)
    ]), "too_many_replicate_lanes")
    monkeypatch.setattr(settings, "MAX_BANDS", 2)
    _expect_422(_payload(lanes=[
        {"name": "lane-1", "observed_bands": [100, 200, 300]},
    ], min_occurrences=1), "too_many_bands")
