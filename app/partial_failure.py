"""核心算法：从异常凝胶反推最可能未切开的酶切边界（部分酶切失败诊断）。

所有函数均为纯函数，不依赖任何外部服务或全局状态。

部分失败合并规则
----------------
请求给出环状质粒按圆环顺序排列的理论片段 f0, f1, …, f_{n-1}（长度 l_i）
及每段的**后边界** b_i（片段 f_i 结束于边界 b_i；末段 f_{n-1} 的后边界
b_{n-1} 跨过坐标原点与 f0 相接）。边界 b_i 若**切开**，f_i 在其末尾断开；
若**失败（未切开）**，f_i 与其后片段在圆环上合并：

  * 连续多个失败边界把其两侧相邻片段并入同一个合并条带；
  * 边界 b_{n-1} 失败时，末段 f_{n-1} 与首段 f0 跨坐标原点合并；
  * 全部边界都失败时，整个圆环合并为 1 条长度等于质粒长度的条带；
  * 没有失败边界时即为完全酶切，预测条带就是各理论片段本身。

枚举与选定
----------
枚举失败边界数 0..``max_failed_boundaries`` 的**全部**边界子集（边界按
请求中的圆环顺序编号，结果与提交顺序无关），对每个子集：

1. 按上述规则合并相邻片段，得到合并后预测条带（按长度升序后参与匹配）；
2. 复用 :func:`app.identify.match_sample` 的一对一最小代价匹配
   （先最大化匹配对数，再最小化总误差），得到逐条配对、缺失/额外条带。

全部评估子集依次按以下次序比较，最后以边界名序列字典序破同分，
因此诊断唯一且确定：

  1. 未匹配条带总数（缺失 + 额外）升序；
  2. 总误差（各配对长度差之和）升序；
  3. 失败边界数升序；
  4. 排序后的失败边界名序列字典序升序。
"""
from itertools import combinations

from app.identify import match_sample
from app.models import (
    PartialFailureMergedBandOut,
    PartialFailureRequest,
    PartialFailureResponse,
    SampleMatchOut,
)


def merge_failed_fragments(
    lengths: list[int],
    boundary_names: list[str],
    failed: frozenset[int],
) -> list[PartialFailureMergedBandOut]:
    """按失败边界集合合并圆环上的相邻理论片段，返回合并条带（按圆环顺序）。

    ``failed`` 为失败边界的下标集合：边界 ``i`` 失败意味着片段 ``i`` 末尾
    不断开、与片段 ``(i+1) mod n`` 合并。结果的首条合并条带必含第 1 个
    理论片段（下标 0），因此条带按圆环顺序唯一定位。
    """
    n = len(lengths)
    bands: list[PartialFailureMergedBandOut] = []

    def close_band(indices: list[int], inside_override: list[str] | None = None) -> None:
        total = sum(lengths[i] for i in indices)
        # 条带内被跨越的失败边界：除条带最后一个片段外，各组成片段的后边界
        inside = (
            [boundary_names[i] for i in indices[:-1]]
            if inside_override is None else inside_override
        )
        bands.append(PartialFailureMergedBandOut(
            length=float(total),
            fragment_indexes=[i + 1 for i in indices],
            fragment_lengths=[lengths[i] for i in indices],
            failed_boundaries_inside=inside,
            crosses_origin=indices[-1] == n - 1,
        ))

    if len(failed) == n:
        # 全部边界失败：整个圆环合并为一条（含跨原点的末段，全部边界均被跨越）
        close_band(list(range(n)), inside_override=list(boundary_names))
        return bands

    indices: list[int] = []
    for i in range(n):
        indices.append(i)
        if i not in failed:
            # 边界 i 切开：当前合并条带在此结束
            close_band(indices)
            indices = []
    # 未闭合的缓冲只可能由“末段后边界失败（跨原点）”产生：
    # 这些片段经原点并入首条合并条带（首条保证含片段 0）。
    if indices:
        first = bands[0]
        first.length = float(first.length + sum(lengths[i] for i in indices))
        first.fragment_indexes = [i + 1 for i in indices] + first.fragment_indexes
        first.fragment_lengths = [lengths[i] for i in indices] + first.fragment_lengths
        # 跨原点并入的片段其各自后边界也在该条带内部被跨越
        first.failed_boundaries_inside = (
            [boundary_names[i] for i in indices] + first.failed_boundaries_inside
        )
        first.crosses_origin = True
    return bands


def resolve_partial_failure(req: PartialFailureRequest) -> PartialFailureResponse:
    """枚举不超过上限的失败边界子集，按固定次序选定唯一诊断。"""
    lengths = [f.length for f in req.fragments]
    boundary_names = [f.following_boundary for f in req.fragments]
    n = len(lengths)
    observed = list(req.observed_bands)
    tolerance = req.band_match_tolerance

    limit = min(req.max_failed_boundaries, n)

    best_key: tuple[int, float, int, tuple[str, ...]] | None = None
    best_failed: frozenset[int] = frozenset()
    best_match: SampleMatchOut | None = None
    best_bands: list[PartialFailureMergedBandOut] = []
    evaluated = 0

    for size in range(0, limit + 1):
        for subset in combinations(range(n), size):
            evaluated += 1
            failed = frozenset(subset)
            bands = merge_failed_fragments(lengths, boundary_names, failed)
            predicted = sorted(b.length for b in bands)
            match = match_sample(predicted, observed, tolerance)
            unmatched = match.missing_count + match.extra_count
            name_tuple = tuple(sorted(boundary_names[i] for i in subset))
            key = (unmatched, match.total_error, size, name_tuple)
            if best_key is None or key < best_key:
                best_key = key
                best_failed = failed
                best_match = match
                best_bands = bands

    unmatched, total_error, size, name_tuple = best_key
    reason = (
        f"在 {evaluated} 个失败边界数不超过 {req.max_failed_boundaries} 的边界子集中，"
        f"选定失败边界 {list(name_tuple)}（{size} 个）：未匹配条带 {unmatched} 条、"
        f"总误差 {total_error:g} bp、失败边界数 {size}，依次按未匹配条带数升序、"
        f"总误差升序、失败边界数升序、排序后的边界名序列字典序比较后为唯一诊断；"
        + ("全部边界均切开（完全酶切）即与观测最吻合。" if size == 0
           else "上述边界未切开，其相邻理论片段在圆环上（含跨原点情形）合并。")
    )
    return PartialFailureResponse(
        plasmid_length=req.plasmid_length,
        band_match_tolerance=tolerance,
        observed_band_count=len(observed),
        fragment_count=n,
        max_failed_boundaries=req.max_failed_boundaries,
        evaluated_subset_count=evaluated,
        failed_boundary_count=size,
        failed_boundaries=list(name_tuple),
        merged_predicted_bands=sorted(b.length for b in best_bands),
        merged_band_details=best_bands,
        match=best_match,
        selection_reason=reason,
    )
