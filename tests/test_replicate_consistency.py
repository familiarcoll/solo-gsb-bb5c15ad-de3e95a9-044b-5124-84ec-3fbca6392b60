"""凝胶重复实验一致性评估：全局聚类、共识谱与稳定性评分的单元 + API 端到端测试。"""
from fastapi.testclient import TestClient

from app.main import app
from app.models import ReplicateConsistencyRequest
from app.replicate_consistency import (
    assess_replicate_consistency,
    cluster_replicate_bands,
)

client = TestClient(app)


# ---------------------------------------------------------------- 全局聚类（纯函数）

def test_cluster_single_linkage_chains_across_lanes():
    # 100(泳道0) -> 108(泳道1) -> 116(泳道0)：逐级差 8 ≤ 10，链式归入同一簇
    pooled = [(100.0, 0), (116.0, 0), (108.0, 1)]
    clusters = cluster_replicate_bands(pooled, tolerance=10)
    assert len(clusters) == 1
    assert [length for length, _ in clusters[0]] == [100.0, 108.0, 116.0]


def test_cluster_tie_breaks_by_length_then_lane_index():
    # 等长条带并列时按泳道序号升序排列，与提交顺序无关
    pooled = [(100.0, 2), (100.0, 0), (105.0, 1), (100.0, 1)]
    clusters = cluster_replicate_bands(pooled, tolerance=10)
    assert clusters == [[(100.0, 0), (100.0, 1), (100.0, 2), (105.0, 1)]]


def test_cluster_splits_when_gap_exceeds_tolerance():
    pooled = [(100.0, 0), (108.0, 1), (130.0, 0)]
    clusters = cluster_replicate_bands(pooled, tolerance=10)
    assert [[length for length, _ in c] for c in clusters] == [[100.0, 108.0], [130.0]]
    # 边界含等号：恰好相差容差时不分裂
    assert len(cluster_replicate_bands([(100.0, 0), (110.0, 1)], tolerance=10)) == 1


# ---------------------------------------------------------------- 一致性评估（纯函数）

def _req(lanes, tolerance=10, min_occurrences=2):
    return ReplicateConsistencyRequest(
        lanes=[{"name": n, "observed_bands": b} for n, b in lanes],
        band_match_tolerance=tolerance,
        min_occurrences=min_occurrences,
    )


def test_assess_basic_consensus_and_missing():
    req = _req([
        ("lane-1", [230, 320, 450]),
        ("lane-2", [232, 318, 450]),
        ("lane-3", [231, 450]),
    ])
    res = assess_replicate_consistency(req)

    assert res.lane_count == 3
    assert res.total_band_count == 8
    assert res.cluster_count == 3
    assert res.stable_cluster_count == 3
    assert res.low_frequency_cluster_count == 0
    # 共识谱：每簇中位长度（偶数个取中间两值均值）
    assert res.consensus_bands == [231.0, 319.0, 450.0]

    c1, c2, c3 = res.clusters
    assert (c1.cluster_id, c1.median_length, c1.band_count, c1.lane_count) == (1, 231.0, 3, 3)
    assert c1.lanes == ["lane-1", "lane-2", "lane-3"]
    assert [(m.lane, m.length) for m in c1.members] == [
        ("lane-1", 230), ("lane-3", 231), ("lane-2", 232),
    ]
    assert not c1.is_low_frequency
    assert (c2.cluster_id, c2.median_length, c2.lane_count) == (2, 319.0, 2)
    assert c2.lanes == ["lane-1", "lane-2"]
    assert (c3.cluster_id, c3.median_length, c3.lane_count) == (3, 450.0, 3)

    # 逐泳道报告按泳道名字典序
    assert [r.lane for r in res.lanes] == ["lane-1", "lane-2", "lane-3"]
    reports = {r.lane: r for r in res.lanes}

    r1 = reports["lane-1"]
    assert r1.band_count == 3 and r1.matched_band_count == 3
    assert r1.missing_count == 0 and r1.extra_count == 0 and r1.low_frequency_count == 0
    assert r1.anomaly_count == 0
    assert r1.consensus_deviation == 2  # |230-231| + |320-319| + 0
    assert r1.stability_score == 84.62  # 100 / (1 + 2/11)
    assert r1.description

    r2 = reports["lane-2"]
    assert r2.consensus_deviation == 2  # |232-231| + |318-319| + 0
    assert r2.stability_score == 84.62

    r3 = reports["lane-3"]
    assert r3.matched_band_count == 2
    assert r3.missing_bands == [319.0]  # 未参与稳定簇 2，取簇中位长度
    assert r3.missing_count == 1
    assert r3.consensus_deviation == 0
    assert r3.anomaly_count == 1
    assert r3.stability_score == 50.0  # 100 / (1 + 1)

    assert res.most_stable_lanes == ["lane-1", "lane-2"]
    assert res.stability_reason


def test_assess_marks_extra_and_low_frequency_bands():
    req = _req([
        ("lane-1", [100, 102, 500]),
        ("lane-2", [101, 500]),
        ("lane-3", [500, 900]),
    ], tolerance=5)
    res = assess_replicate_consistency(req)

    assert res.cluster_count == 3
    assert res.stable_cluster_count == 2
    assert res.low_frequency_cluster_count == 1
    assert res.consensus_bands == [101.0, 500.0, 900.0]

    c1, _, c3 = res.clusters
    # 100/101/102 链式同簇：3 条带、2 条泳道
    assert (c1.median_length, c1.band_count, c1.lane_count) == (101.0, 3, 2)
    assert not c1.is_low_frequency
    assert c3.is_low_frequency and c3.lane_count == 1

    reports = {r.lane: r for r in res.lanes}
    r1 = reports["lane-1"]
    # 100 与 102 距中位 101 同为 1，并列取较短者 100 为代表，102 计为额外
    assert r1.matched_band_count == 2
    assert r1.extra_bands == [102]
    assert r1.consensus_deviation == 1  # |100-101| + |500-500|
    assert r1.anomaly_count == 1
    assert r1.stability_score == 46.15  # 100 / (1 + 1 + 1/6)

    r2 = reports["lane-2"]
    assert r2.anomaly_count == 0
    assert r2.consensus_deviation == 0
    assert r2.stability_score == 100.0

    r3 = reports["lane-3"]
    # 未参与稳定簇 [100,101,102] 计缺失（取中位 101），900 落入低频簇计低频
    assert r3.missing_bands == [101.0]
    assert r3.missing_count == 1
    assert r3.low_frequency_bands == [900]
    assert r3.low_frequency_count == 1
    assert r3.anomaly_count == 2
    assert r3.stability_score == 33.33  # 100 / (1 + 2)

    # 低频簇不产生缺失：lane-1/lane-2 未覆盖 900 簇但不计缺失
    assert r1.missing_count == 0 and r2.missing_count == 0
    assert res.most_stable_lanes == ["lane-2"]


def test_assess_min_occurrences_one_marks_every_cluster_stable():
    req = _req([
        ("lane-1", [100, 500]),
        ("lane-2", [500, 900]),
    ], tolerance=5, min_occurrences=1)
    res = assess_replicate_consistency(req)
    assert res.low_frequency_cluster_count == 0
    assert res.stable_cluster_count == 3
    # 每条泳道都缺失对方独有的稳定簇
    reports = {r.lane: r for r in res.lanes}
    assert reports["lane-1"].missing_bands == [900.0]
    assert reports["lane-2"].missing_bands == [100.0]


def test_assess_single_lane_all_bands_stable():
    req = _req([("only", [100, 102, 500])], tolerance=5, min_occurrences=1)
    res = assess_replicate_consistency(req)
    assert res.cluster_count == 2
    report = res.lanes[0]
    assert report.missing_count == 0 and report.low_frequency_count == 0
    # 100/102 同簇，中位 101，代表 100，102 计为额外
    assert report.extra_bands == [102]
    assert report.consensus_deviation == 1
    assert res.most_stable_lanes == ["only"]


def test_assess_is_order_independent():
    first = assess_replicate_consistency(_req([
        ("lane-b", [232, 318, 450]),
        ("lane-a", [230, 320, 450]),
        ("lane-c", [231, 450]),
    ]))
    second = assess_replicate_consistency(_req([
        ("lane-c", [450, 231]),
        ("lane-a", [450, 320, 230]),
        ("lane-b", [318, 450, 232]),
    ]))
    assert first == second


# ---------------------------------------------------------------- API 端到端

def _payload(**kw):
    payload = {
        "lanes": [
            {"name": "lane-1", "observed_bands": [230, 320, 450]},
            {"name": "lane-2", "observed_bands": [232, 318, 450]},
            {"name": "lane-3", "observed_bands": [231, 450]},
        ],
        "band_match_tolerance": 10,
        "min_occurrences": 2,
    }
    payload.update(kw)
    return payload


def test_api_replicate_consistency_basic():
    res = client.post("/api/v1/digest/replicate-consistency", json=_payload())
    assert res.status_code == 200
    body = res.json()

    assert body["band_match_tolerance"] == 10
    assert body["min_occurrences"] == 2
    assert body["lane_count"] == 3
    assert body["total_band_count"] == 8
    assert body["cluster_count"] == 3
    assert body["stable_cluster_count"] == 3
    assert body["low_frequency_cluster_count"] == 0
    assert body["consensus_bands"] == [231.0, 319.0, 450.0]
    assert body["most_stable_lanes"] == ["lane-1", "lane-2"]
    assert body["stability_reason"]

    # 簇组成依据：编号从 1 开始、按中位长度升序
    assert [c["cluster_id"] for c in body["clusters"]] == [1, 2, 3]
    for c in body["clusters"]:
        assert c["band_count"] == len(c["members"])
        assert c["lane_count"] == len(c["lanes"])
        assert c["min_length"] <= c["median_length"] <= c["max_length"]
        # members 按（长度、泳道序号）升序
        keys = [(m["length"], m["lane"]) for m in c["members"]]
        assert keys == sorted(keys)

    # 逐泳道计数自洽：条带数 = 代表数 + 额外数；异常数 = 缺失 + 额外 + 低频
    for report in body["lanes"]:
        assert report["band_count"] == report["matched_band_count"] + report["extra_count"]
        assert report["anomaly_count"] == (
            report["missing_count"] + report["extra_count"] + report["low_frequency_count"]
        )
        assert report["missing_count"] == len(report["missing_bands"])
        assert report["extra_count"] == len(report["extra_bands"])
        assert report["low_frequency_count"] == len(report["low_frequency_bands"])
        assert report["description"]


def test_api_replicate_consistency_order_independent():
    first = client.post("/api/v1/digest/replicate-consistency", json=_payload()).json()
    shuffled = _payload(lanes=[
        {"name": "lane-3", "observed_bands": [450, 231]},
        {"name": "lane-1", "observed_bands": [450, 320, 230]},
        {"name": "lane-2", "observed_bands": [318, 450, 232]},
    ])
    second = client.post("/api/v1/digest/replicate-consistency", json=shuffled).json()
    assert first == second


def test_api_replicate_consistency_float_bands():
    body = client.post("/api/v1/digest/replicate-consistency", json=_payload(
        lanes=[
            {"name": "lane-1", "observed_bands": [100.4, 500]},
            {"name": "lane-2", "observed_bands": [99.6, 500]},
        ],
        band_match_tolerance=1,
    )).json()
    assert body["consensus_bands"] == [100.0, 500.0]
    for report in body["lanes"]:
        assert report["anomaly_count"] == 0
        assert report["consensus_deviation"] == 0.4


# ---------------------------------------------------------------- 422 校验

def _expect_replicate_422(payload, code_part=None):
    res = client.post("/api/v1/digest/replicate-consistency", json=payload)
    assert res.status_code == 422
    body = res.json()
    assert body["code"] == "validation_failed"
    assert body["errors"]
    if code_part:
        assert any(code_part in e["type"] for e in body["errors"])
    return body


def test_replicate_422_empty_lanes_and_empty_bands():
    body = _expect_replicate_422(_payload(lanes=[]))
    assert any("lanes" in e["loc"] for e in body["errors"])
    _expect_replicate_422(_payload(lanes=[{"name": "lane-1", "observed_bands": []}]))


def test_replicate_422_blank_and_duplicate_lane_names():
    _expect_replicate_422(_payload(lanes=[
        {"name": "   ", "observed_bands": [100]},
        {"name": "lane-2", "observed_bands": [100]},
    ]), "invalid_lane_name")
    _expect_replicate_422(_payload(lanes=[
        {"name": "lane-1", "observed_bands": [100]},
        {"name": "lane-1", "observed_bands": [200]},
    ]), "duplicate_lane_name")


def test_replicate_422_negative_band_and_bad_tolerance():
    _expect_replicate_422(_payload(lanes=[{"name": "lane-1", "observed_bands": [-1]}]))
    _expect_replicate_422(_payload(band_match_tolerance=-1))
    _expect_replicate_422(_payload(band_match_tolerance="wide"))
    _expect_replicate_422(_payload(band_match_tolerance=None))


def test_replicate_422_min_occurrences_invalid():
    _expect_replicate_422(_payload(min_occurrences=0))
    _expect_replicate_422(_payload(min_occurrences=-1))
    _expect_replicate_422(_payload(min_occurrences="two"))
    body = _expect_replicate_422(_payload(min_occurrences=4), "min_occurrences_too_large")
    assert any("min_occurrences" in e["loc"] for e in body["errors"])


def test_replicate_422_missing_required_fields():
    # 缺泳道列表、缺容差、缺最少出现次数、泳道缺名称或缺条带
    _expect_replicate_422({})
    _expect_replicate_422({"band_match_tolerance": 10, "min_occurrences": 1})
    _expect_replicate_422(_payload(lanes=[{"observed_bands": [100]}]))
    _expect_replicate_422(_payload(lanes=[{"name": "lane-1"}]))


def test_replicate_422_too_many_lanes(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_REPLICATE_LANES", 2)
    _expect_replicate_422(_payload(), "too_many_replicate_lanes")


def test_replicate_422_too_many_bands(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_BANDS", 2)
    _expect_replicate_422(_payload(), "too_many_bands")


def test_replicate_422_reports_multiple_errors_at_once():
    # 重名 + 最少出现次数超泳道数：一次性报告
    body = _expect_replicate_422(_payload(
        lanes=[
            {"name": "lane-1", "observed_bands": [100]},
            {"name": "lane-1", "observed_bands": [200]},
        ],
        min_occurrences=5,
    ))
    types = {e["type"] for e in body["errors"]}
    assert "duplicate_lane_name" in types
    assert "min_occurrences_too_large" in types
