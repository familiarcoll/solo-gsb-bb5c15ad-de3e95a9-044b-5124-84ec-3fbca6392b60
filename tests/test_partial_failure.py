"""部分酶切失败反推：失败边界合并、子集枚举与唯一诊断选定的单元 + API 端到端测试。"""
from fastapi.testclient import TestClient

from app.main import app
from app.models import PartialFailureRequest
from app.partial_failure import merge_failed_fragments, resolve_partial_failure

client = TestClient(app)


# ---------------------------------------------------------------- 失败边界合并（纯函数）

LENGTHS = [200, 300, 100, 400]
NAMES = ["B1", "B2", "B3", "B4"]


def _bands(failed):
    return merge_failed_fragments(LENGTHS, NAMES, frozenset(failed))


def test_merge_no_failure_is_identity():
    bands = _bands([])
    assert [b.length for b in bands] == [200.0, 300.0, 100.0, 400.0]
    assert [b.fragment_indexes for b in bands] == [[1], [2], [3], [4]]
    assert sum(b.length for b in bands) == 1000.0
    # 仅末段条带跨原点
    assert [b.crosses_origin for b in bands] == [False, False, False, True]
    assert all(b.failed_boundaries_inside == [] for b in bands)


def test_merge_single_failed_boundary_joins_adjacent_fragments():
    # B2 失败：f1 与 f2 合并为 400
    bands = _bands([1])
    assert [b.length for b in bands] == [200.0, 400.0, 400.0]
    merged = bands[1]
    assert merged.fragment_indexes == [2, 3]
    assert merged.fragment_lengths == [300, 100]
    assert merged.failed_boundaries_inside == ["B2"]
    assert merged.crosses_origin is False


def test_merge_consecutive_failures_form_one_chain():
    # B1、B2 连续失败：f0+f1+f2 合并为 600
    bands = _bands([0, 1])
    assert [b.length for b in bands] == [600.0, 400.0]
    assert bands[0].fragment_indexes == [1, 2, 3]
    assert bands[0].failed_boundaries_inside == ["B1", "B2"]
    assert bands[0].crosses_origin is False
    assert bands[1].crosses_origin is True


def test_merge_last_boundary_failure_crosses_origin():
    # B4（末段后边界）失败：f3 与 f0 跨原点合并为 600，首条即跨原点条带
    bands = _bands([3])
    first = bands[0]
    assert first.length == 600.0
    assert first.fragment_indexes == [4, 1]
    assert first.fragment_lengths == [400, 200]
    assert first.failed_boundaries_inside == ["B4"]
    assert first.crosses_origin is True
    assert [b.length for b in bands] == [600.0, 300.0, 100.0]


def test_merge_consecutive_failures_across_origin():
    # B4 与 B1 失败在圆环上连续（跨原点）：f3+f0+f1 合并为 900
    bands = _bands([3, 0])
    assert [b.length for b in bands] == [900.0, 100.0]
    first = bands[0]
    assert first.fragment_indexes == [4, 1, 2]
    assert first.fragment_lengths == [400, 200, 300]
    assert first.failed_boundaries_inside == ["B4", "B1"]
    assert first.crosses_origin is True


def test_merge_all_boundaries_fail_gives_single_full_length_band():
    bands = _bands([0, 1, 2, 3])
    assert len(bands) == 1
    only = bands[0]
    assert only.length == 1000.0
    assert only.fragment_indexes == [1, 2, 3, 4]
    assert only.failed_boundaries_inside == ["B1", "B2", "B3", "B4"]
    assert only.crosses_origin is True


def test_merge_single_fragment_circular_plasmid():
    # n = 1：唯一边界失败即整环（长度 = 质粒长度），切开则为该片段本身
    cut = merge_failed_fragments([500], ["B0"], frozenset())
    assert [b.length for b in cut] == [500.0]
    assert cut[0].crosses_origin is True
    failed = merge_failed_fragments([500], ["B0"], frozenset([0]))
    assert [b.length for b in failed] == [500.0]
    assert failed[0].failed_boundaries_inside == ["B0"]


def test_merge_band_lengths_always_sum_to_plasmid_length():
    for failed in ([], [0], [1], [2], [3], [0, 1], [2, 3], [3, 0], [0, 2],
                   [0, 1, 2], [0, 1, 3], [0, 1, 2, 3]):
        bands = _bands(failed)
        assert sum(b.length for b in bands) == 1000.0, failed
        indexes = sorted(i for b in bands for i in b.fragment_indexes)
        assert indexes == [1, 2, 3, 4], failed


# ---------------------------------------------------------------- 诊断选定（纯函数）

def _req(observed, tolerance=10, max_failed=2, lengths=None, names=None):
    lengths = LENGTHS if lengths is None else lengths
    names = NAMES if names is None else names
    return PartialFailureRequest(
        plasmid_length=sum(lengths),
        fragments=[
            {"length": l, "following_boundary": n}
            for l, n in zip(lengths, names)
        ],
        observed_bands=observed,
        band_match_tolerance=tolerance,
        max_failed_boundaries=max_failed,
    )


def test_resolve_diagnoses_single_failed_boundary():
    # 完全酶切预测为 [200,300,100,400]；B2 失败合并出 [200,400,400]，
    # 观测恰为 3 条带 [201,398,400] -> B2 失败为唯一零未匹配诊断
    res = resolve_partial_failure(_req([201, 398, 400]))
    assert res.failed_boundaries == ["B2"]
    assert res.failed_boundary_count == 1
    assert res.merged_predicted_bands == [200.0, 400.0, 400.0]
    assert res.evaluated_subset_count == 1 + 4 + 6  # C(4,0)+C(4,1)+C(4,2)
    match = res.match
    assert match.matched_count == 3
    assert match.missing_count == 0 and match.extra_count == 0
    assert match.total_error == 3  # 1 + 2 + 0
    # 合并明细与预测条带自洽
    assert sorted(b.length for b in res.merged_band_details) == res.merged_predicted_bands


def test_resolve_complete_digest_wins_when_observation_matches():
    res = resolve_partial_failure(_req([200, 300, 100, 400]))
    assert res.failed_boundaries == []
    assert res.failed_boundary_count == 0
    assert res.match.total_error == 0
    assert res.merged_predicted_bands == [100.0, 200.0, 300.0, 400.0]


def test_resolve_prefers_fewer_failed_boundaries_on_tie():
    # 对称四段（各 100），观测 [100, 200]、容差 0：
    #   0 失败预测 [100,100,100,100] -> 匹配 1（100），未匹配 3
    #   单失败预测 [100,100,200]     -> 匹配 2，未匹配 1
    #   非相邻两失败（如 {A,C}）预测 [200,200] -> 仅 1 条 200 可配对，匹配 1，未匹配 1
    # 单失败与两失败未匹配数相同（各 1 条），按失败边界数更少取单失败
    lengths = [100, 100, 100, 100]
    req = _req([100, 200], tolerance=0, max_failed=2, lengths=lengths,
               names=["A", "B", "C", "D"])
    res = resolve_partial_failure(req)
    assert res.failed_boundary_count == 1
    assert res.match.missing_count + res.match.extra_count == 1


def test_resolve_total_error_tiebreak():
    # 单边界失败中只有 B1 能把预测合并为 3 条带解释 3 条观测：
    # B1 失败 -> [500,100,400]，观测 [502,99,402]，总误差 5
    req = _req([502, 99, 402], tolerance=10)
    res = resolve_partial_failure(req)
    assert res.failed_boundaries == ["B1"]
    assert res.match.total_error == 5
    assert res.match.matched_count == 3


def test_resolve_final_tiebreak_is_boundary_name_sequence():
    # 对称圆环（四段等长）：任何单边界失败的预测完全相同，
    # 最终按排序后的边界名序列字典序取第一个
    lengths = [250, 250, 250, 250]
    req = _req([500, 250, 250], tolerance=0, max_failed=1,
               lengths=lengths, names=["zeta", "alpha", "mid", "beta"])
    res = resolve_partial_failure(req)
    assert res.failed_boundary_count == 1
    assert res.failed_boundaries == ["alpha"]
    assert res.selection_reason


def test_resolve_max_failed_boundaries_limits_enumeration():
    req = _req([600, 400], max_failed=1)
    # 只枚举 0..1 个失败：5 个子集；连续两失败（600+400）不参与
    assert resolve_partial_failure(req).evaluated_subset_count == 5
    full = resolve_partial_failure(req.model_copy(update={"max_failed_boundaries": 2}))
    assert full.evaluated_subset_count == 11
    assert full.failed_boundaries == ["B1", "B2"]
    assert full.merged_predicted_bands == [400.0, 600.0]
    assert full.match.total_error == 0


def test_resolve_all_boundaries_fail_when_limit_covers_all():
    # 观测只有一条全长带；3 个连续失败已能把整环合并为 1000（4 个失败同样
    # 只剩一条 1000 带，指标相同），按失败边界数更少取 3 个
    req = _req([1000], max_failed=4)
    res = resolve_partial_failure(req)
    assert res.failed_boundary_count == 3
    assert res.merged_predicted_bands == [1000.0]
    assert res.match.matched_count == 1
    # 上限允许时 4 个边界全失败子集同样被评估（共 16 个子集）
    assert res.evaluated_subset_count == 16


def test_resolve_orders_failed_boundaries_by_name_but_details_keep_ring_order():
    req = _req([900, 100], max_failed=2)
    res = resolve_partial_failure(req)
    # B4、B1 跨原点连续失败合并出 900
    assert res.failed_boundaries == ["B1", "B4"]  # 名称字典序
    first = res.merged_band_details[0]
    assert first.failed_boundaries_inside == ["B4", "B1"]  # 圆环顺序
    assert first.crosses_origin is True


# ---------------------------------------------------------------- API 端到端

def _payload(**kw):
    payload = {
        "plasmid_length": 1000,
        "fragments": [
            {"length": 200, "following_boundary": "B1"},
            {"length": 300, "following_boundary": "B2"},
            {"length": 100, "following_boundary": "B3"},
            {"length": 400, "following_boundary": "B4"},
        ],
        "observed_bands": [201, 398, 400],
        "band_match_tolerance": 10,
        "max_failed_boundaries": 2,
    }
    payload.update(kw)
    return payload


def test_api_partial_failure_full_evidence():
    res = client.post("/api/v1/digest/partial-failure", json=_payload())
    assert res.status_code == 200
    body = res.json()

    assert body["fragment_count"] == 4
    assert body["max_failed_boundaries"] == 2
    assert body["evaluated_subset_count"] == 11
    assert body["failed_boundaries"] == ["B2"]
    assert body["failed_boundary_count"] == 1
    assert body["observed_band_count"] == 3
    assert body["band_match_tolerance"] == 10
    assert body["merged_predicted_bands"] == [200.0, 400.0, 400.0]
    assert body["selection_reason"]

    match = body["match"]
    assert match["matched_count"] == 3
    assert match["missing_count"] == 0 and match["extra_count"] == 0
    assert match["missing_bands"] == [] and match["extra_bands"] == []
    assert match["total_error"] == 3
    assert match["predicted_band_count"] == len(body["merged_predicted_bands"])
    for m in match["matches"]:
        assert m["error"] == abs(m["predicted_length"] - m["observed_length"])
        assert m["description"]

    # 合并明细：组成片段长度之和 = 合并条带长度，且总长守恒
    for band in body["merged_band_details"]:
        assert band["length"] == sum(band["fragment_lengths"])
    assert sum(b["length"] for b in body["merged_band_details"]) == 1000
    # 恰好一个条带包含被跨越的失败边界 B2
    failed_inside = [
        b for b in body["merged_band_details"] if b["failed_boundaries_inside"]
    ]
    assert len(failed_inside) == 1
    assert failed_inside[0]["failed_boundaries_inside"] == ["B2"]


def test_api_partial_failure_cross_origin():
    body = client.post("/api/v1/digest/partial-failure", json=_payload(
        observed_bands=[600, 300, 100],
    )).json()
    assert body["failed_boundaries"] == ["B4"]
    first = body["merged_band_details"][0]
    assert first["crosses_origin"] is True
    assert first["fragment_indexes"] == [4, 1]
    assert first["failed_boundaries_inside"] == ["B4"]


def test_api_partial_failure_counts_consistent():
    body = client.post("/api/v1/digest/partial-failure", json=_payload(
        observed_bands=[201, 999],
    )).json()
    match = body["match"]
    assert match["matched_count"] + match["missing_count"] == match["predicted_band_count"]
    assert match["matched_count"] + match["extra_count"] == body["observed_band_count"]
    assert match["missing_count"] == len(match["missing_bands"])
    assert match["extra_count"] == len(match["extra_bands"])


# ---------------------------------------------------------------- 422 校验

def _expect_422(payload, code_part=None):
    res = client.post("/api/v1/digest/partial-failure", json=payload)
    assert res.status_code == 422
    body = res.json()
    assert body["code"] == "validation_failed"
    assert body["errors"]
    if code_part:
        assert any(code_part in e["type"] for e in body["errors"])
    return body


def test_422_fragment_total_mismatch():
    body = _expect_422(_payload(plasmid_length=999), "fragment_length_mismatch")
    assert any("fragments" in e["loc"] for e in body["errors"])


def test_422_duplicate_boundary_name():
    fragments = [
        {"length": 200, "following_boundary": "B1"},
        {"length": 300, "following_boundary": "B1"},
        {"length": 100, "following_boundary": "B3"},
        {"length": 400, "following_boundary": "B4"},
    ]
    _expect_422(_payload(fragments=fragments), "duplicate_boundary_name")


def test_422_blank_boundary_name():
    fragments = [
        {"length": 200, "following_boundary": "  "},
        {"length": 800, "following_boundary": "B2"},
    ]
    _expect_422(_payload(fragments=fragments), "invalid_boundary_name")


def test_422_nonpositive_fragment_length():
    fragments = [
        {"length": 0, "following_boundary": "B1"},
        {"length": 1000, "following_boundary": "B2"},
    ]
    _expect_422(_payload(fragments=fragments))


def test_422_max_failed_boundaries_above_fragment_count():
    body = _expect_422(_payload(max_failed_boundaries=5),
                       "max_failed_boundaries_too_large")
    assert any("max_failed_boundaries" in e["loc"] for e in body["errors"])


def test_422_max_failed_boundaries_negative_or_missing():
    _expect_422(_payload(max_failed_boundaries=-1))
    payload = _payload()
    del payload["max_failed_boundaries"]
    _expect_422(payload)


def test_422_empty_fragments_and_empty_observed():
    _expect_422(_payload(fragments=[]))
    _expect_422(_payload(observed_bands=[]))


def test_422_negative_observed_band_and_bad_tolerance():
    _expect_422(_payload(observed_bands=[100, -5]))
    _expect_422(_payload(band_match_tolerance=-1))
    _expect_422(_payload(band_match_tolerance=None))


def test_422_too_many_boundaries(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_PARTIAL_FAILURE_BOUNDARIES", 3)
    fragments = [
        {"length": 250, "following_boundary": f"B{i}"} for i in range(4)
    ]
    _expect_422(_payload(fragments=fragments), "too_many_partial_failure_boundaries")


def test_422_reports_multiple_errors_at_once():
    # 总长不符 + max_failed_boundaries 超过边界数 + 边界重名，一次性报告
    fragments = [
        {"length": 200, "following_boundary": "B1"},
        {"length": 300, "following_boundary": "B1"},
        {"length": 100, "following_boundary": "B3"},
        {"length": 400, "following_boundary": "B4"},
    ]
    body = _expect_422(_payload(
        plasmid_length=998, fragments=fragments, max_failed_boundaries=9,
    ))
    types = {e["type"] for e in body["errors"]}
    assert "fragment_length_mismatch" in types
    assert "duplicate_boundary_name" in types
    assert "max_failed_boundaries_too_large" in types
