"""核心算法：同一样本多泳道重复实验的全局条带聚类与一致性评估。

所有函数均为纯函数，不依赖任何外部服务或全局状态。

全局条带聚类
------------
请求给出同一样本的若干**命名观测泳道**、条带匹配容差与最少出现次数
``min_occurrences``。先把泳道按名称字典序编号（泳道序号），再把全部泳道的
观测条带汇总后按 **(长度, 泳道序号)** 升序排序——并列按长度、再按泳道序号
破同分——做单连接聚类：相邻条带长度差 ≤ 容差即归入同一簇（链式可传递）。
排序键只依赖条带长度与泳道名，因此聚类结果**与输入提交顺序无关**。

共识谱与逐泳道标记
------------------
每个簇以其全部成员条带长度的**中位数**（偶数个取中间两值均值，保留两位小数）
作为代表长度，全部簇的中位长度按升序构成**共识谱**。簇内不同泳道数为该条带的
**出现次数**：达到 ``min_occurrences`` 的簇为**稳定簇**，否则为**低频簇**。

对每条泳道：

* **缺失条带**：该泳道未参与的稳定簇，取其簇中位长度；
* **额外条带**：该泳道在同一簇内有多条条带时，距簇中位长度最近
  （并列取较短者）的一条作为代表条带，其余计为额外条带；
* **低频条带**：该泳道的代表条带落入低频簇。

**共识偏差**为该泳道各簇代表条带与对应簇中位长度差绝对值之和（bp）；
**异常数** = 缺失 + 额外 + 低频。

稳定性评分
----------
``stability_score = 100 / (1 + 异常数 + 共识偏差 / (容差 + 1))``，保留两位小数，
满分 100（无异常且零偏差）；共识偏差以容差为尺度折算后与异常数同量纲相加，
评分只使用响应中输出的计数值与（保留两位小数的）共识偏差，可由响应直接复算。
"""
from app.models import (
    ReplicateClusterMemberOut,
    ReplicateClusterOut,
    ReplicateConsistencyRequest,
    ReplicateConsistencyResponse,
    ReplicateLaneReportOut,
)


def _median_length(lengths: list[float]) -> float:
    """簇内全部条带长度的中位数：偶数个取中间两值的均值，保留两位小数。"""
    ordered = sorted(lengths)
    n = len(ordered)
    mid = n // 2
    if n % 2 == 1:
        median = ordered[mid]
    else:
        median = (ordered[mid - 1] + ordered[mid]) / 2
    return round(median, 2)


def cluster_replicate_bands(pooled: list[tuple[float, int]],
                            tolerance: float) -> list[list[tuple[float, int]]]:
    """把 (条带长度, 泳道序号) 列表排序后做单连接全局聚类。

    排序键为 (长度, 泳道序号)：并列按长度、再按泳道序号破同分；相邻条带
    长度差 ≤ 容差即同簇（链式可传递）。结果与输入提交顺序无关。
    """
    ordered = sorted(pooled, key=lambda item: (item[0], item[1]))
    clusters: list[list[tuple[float, int]]] = []
    for item in ordered:
        if clusters and item[0] - clusters[-1][-1][0] <= tolerance:
            clusters[-1].append(item)
        else:
            clusters.append([item])
    return clusters


def assess_replicate_consistency(
    req: ReplicateConsistencyRequest,
) -> ReplicateConsistencyResponse:
    """全局聚类所有泳道的观测条带，形成共识谱并逐泳道评估一致性。"""
    lane_names = sorted(w.name for w in req.lanes)
    bands_by_lane = {w.name: list(w.observed_bands) for w in req.lanes}
    lane_index = {name: i for i, name in enumerate(lane_names)}
    tolerance = req.band_match_tolerance
    min_occurrences = req.min_occurrences

    # 1) 全局单连接聚类：泳道按名字典序编号，并列按长度、泳道序号破同分
    pooled: list[tuple[float, int]] = []
    for name in lane_names:
        for band in bands_by_lane[name]:
            pooled.append((band, lane_index[name]))
    raw_clusters = cluster_replicate_bands(pooled, tolerance)

    # 2) 每簇统计：中位长度、出现次数（簇内不同泳道数）与低频标记
    clusters: list[ReplicateClusterOut] = []
    # (中位长度, 是否低频, {泳道序号: [该泳道在簇内的条带]})
    cluster_infos: list[tuple[float, bool, dict[int, list[float]]]] = []
    for cid, members in enumerate(raw_clusters, start=1):
        lengths = [length for length, _ in members]
        median = _median_length(lengths)
        by_lane: dict[int, list[float]] = {}
        for length, idx in members:
            by_lane.setdefault(idx, []).append(length)
        lane_count = len(by_lane)
        is_low = lane_count < min_occurrences
        clusters.append(ReplicateClusterOut(
            cluster_id=cid,
            median_length=median,
            min_length=min(lengths),
            max_length=max(lengths),
            band_count=len(lengths),
            lane_count=lane_count,
            is_low_frequency=is_low,
            lanes=[lane_names[i] for i in sorted(by_lane)],
            members=[
                ReplicateClusterMemberOut(lane=lane_names[idx], length=length)
                for length, idx in members
            ],
        ))
        cluster_infos.append((median, is_low, by_lane))

    # 3) 逐泳道标记：代表条带、缺失、额外、低频与共识偏差
    reports: list[ReplicateLaneReportOut] = []
    for idx, name in enumerate(lane_names):
        missing: list[float] = []
        extra: list[float] = []
        low_frequency: list[float] = []
        deviation = 0.0
        matched = 0
        for median, is_low, by_lane in cluster_infos:
            own = by_lane.get(idx)
            if own is None:
                if not is_low:
                    missing.append(median)
                continue
            matched += 1
            # 代表条带：距簇中位长度最近，并列取较短者；其余计为额外条带
            representative = min(own, key=lambda b: (abs(b - median), b))
            rest = sorted(own)
            rest.remove(representative)
            extra.extend(rest)
            deviation += abs(representative - median)
            if is_low:
                low_frequency.append(representative)
        deviation = round(deviation, 2)
        anomaly = len(missing) + len(extra) + len(low_frequency)
        score = round(100 / (1 + anomaly + deviation / (tolerance + 1)), 2)
        reports.append(ReplicateLaneReportOut(
            lane=name,
            band_count=len(bands_by_lane[name]),
            matched_band_count=matched,
            missing_count=len(missing),
            missing_bands=sorted(missing),
            extra_count=len(extra),
            extra_bands=sorted(extra),
            low_frequency_count=len(low_frequency),
            low_frequency_bands=sorted(low_frequency),
            anomaly_count=anomaly,
            consensus_deviation=deviation,
            stability_score=score,
            description=(
                f"泳道 {name}：{len(bands_by_lane[name])} 条观测条带，"
                f"作为代表参与 {matched} 个簇；缺失 {len(missing)} 条、"
                f"额外 {len(extra)} 条、低频 {len(low_frequency)} 条，"
                f"异常数 {anomaly}；共识偏差 {deviation:g} bp，"
                f"稳定性评分 {score}"
            ),
        ))

    best_score = max(r.stability_score for r in reports)
    most_stable = [r.lane for r in reports if r.stability_score == best_score]

    stable_count = sum(1 for c in clusters if not c.is_low_frequency)
    low_count = len(clusters) - stable_count
    total_bands = sum(len(b) for b in bands_by_lane.values())
    reason = (
        f"全部 {len(lane_names)} 条泳道共 {total_bands} 条观测条带，"
        f"按（长度、泳道序号）升序排序后以容差 {tolerance:g} bp 做单连接全局聚类，"
        f"并列按长度与泳道序号破同分，共形成 {len(clusters)} 个簇"
        f"（稳定簇 {stable_count} 个、低频簇 {low_count} 个）；"
        f"以每簇中位长度形成共识谱，出现次数（簇内不同泳道数）达到最少出现次数 "
        f"{min_occurrences} 的簇为稳定簇；逐泳道标记缺失（未覆盖稳定簇）、"
        f"额外（同簇多条带中的非代表条带）与低频（代表条带落入低频簇）条带，"
        f"异常数 = 缺失 + 额外 + 低频，"
        f"稳定性评分 = 100 / (1 + 异常数 + 共识偏差 / (容差 + 1))；"
        f"聚类与评分只依赖条带长度与泳道名，结果与输入顺序无关"
    )

    return ReplicateConsistencyResponse(
        band_match_tolerance=tolerance,
        min_occurrences=min_occurrences,
        lane_count=len(lane_names),
        total_band_count=total_bands,
        cluster_count=len(clusters),
        stable_cluster_count=stable_count,
        low_frequency_cluster_count=low_count,
        consensus_bands=[c.median_length for c in clusters],
        clusters=clusters,
        lanes=reports,
        most_stable_lanes=most_stable,
        stability_reason=reason,
    )
