"""核心算法：凝胶重复实验（同一样本多条命名观测泳道）的一致性评估。

所有函数均为纯函数，不依赖任何外部服务或全局状态。

全局条带聚类
------------
把所有泳道的观测条带汇入全局条带池，每条条带记录 ``(长度, 泳道序号,
道内序号)``：泳道序号为泳道名按字典序排序后的名次，道内序号为该泳道内
条带按长度升序的名次。条带池按 ``(长度, 泳道序号, 道内序号)`` 升序排列
——**并列按长度与泳道序号破同分**，因此聚类结果与泳道、条带的输入提交
顺序无关。扫描有序条带池，相邻条带长度差 ≤ 匹配容差即归入同一簇（等价
于以容差为阈值的单连接聚类，链式可传递）。

共识谱与逐泳道标记
------------------
每个簇的**出现次数**为贡献该簇的不同泳道数（同一泳道的多条条带只计
一次）。出现次数 ≥ ``min_occurrences`` 的簇进入**共识谱**，共识条带长度
取簇内全部条带长度的**中位数**（偶数条取中间两条的均值，保留两位小数）；
其余簇为**低频簇**。对每条泳道：

  * **缺失条带**：该泳道未贡献条带的共识簇（以共识中位长度记录）；
  * **额外条带**：该泳道在同一共识簇内除代表条带（与中位长度最接近、
    并列取较短者）以外的条带；
  * **低频条带**：该泳道落入低频簇的条带。

稳定性评分
----------
每条泳道的共识偏差为其各命中共识簇 ``|代表条带 − 中位长度|`` 之和，
平均共识偏差为共识偏差 ÷ 命中共识条带数（未命中时为 0）。异常数 =
缺失 + 额外 + 低频条带数。稳定性评分（0~100）：

    100 × 命中共识条带数 ÷（共识条带数 + 额外条带数 + 低频条带数）
        ÷（1 + 平均共识偏差）

无异常且共识偏差为 0 时得 100；异常越多、偏差越大，评分越趋近 0。评分
只使用响应中输出的计数值与保留两位小数的偏差，可由响应直接复算；整体
稳定性为各泳道评分（保留两位小数）的均值。
"""
from app.models import (
    ReplicateClusterMemberOut,
    ReplicateClusterOut,
    ReplicateConsistencyRequest,
    ReplicateConsistencyResponse,
    ReplicateLaneConsistencyOut,
)


def _median(sorted_lengths: list[float]) -> float:
    """有序长度列表的中位数（偶数条取中间两条的均值），保留两位小数。"""
    n = len(sorted_lengths)
    mid = n // 2
    if n % 2 == 1:
        return round(sorted_lengths[mid], 2)
    return round((sorted_lengths[mid - 1] + sorted_lengths[mid]) / 2, 2)


def assess_replicate_consistency(
    req: ReplicateConsistencyRequest,
) -> ReplicateConsistencyResponse:
    """全局条带聚类 → 共识谱 → 逐泳道缺失/额外/低频标记与稳定性评分。"""
    lane_names = sorted(lane.name for lane in req.lanes)
    lane_rank = {name: i for i, name in enumerate(lane_names)}
    bands_by_lane = {lane.name: sorted(lane.observed_bands) for lane in req.lanes}
    tolerance = req.band_match_tolerance

    # 1) 全局条带池：(长度, 泳道序号, 道内序号)；并列按长度与泳道序号破同分
    pool: list[tuple[float, int, int]] = []
    for name in lane_names:
        for band_index, length in enumerate(bands_by_lane[name]):
            pool.append((length, lane_rank[name], band_index))
    pool.sort()

    # 2) 确定性的全局单连接聚类：相邻条带长度差 ≤ 容差即同簇（链式可传递）
    clusters: list[list[tuple[float, int, int]]] = []
    for item in pool:
        if clusters and item[0] - clusters[-1][-1][0] <= tolerance:
            clusters[-1].append(item)
        else:
            clusters.append([item])

    # 3) 逐簇统计：出现次数（不同泳道数）、中位长度、是否入共识谱，
    #    并累计各泳道的缺失 / 额外 / 低频条带与共识偏差
    cluster_outs: list[ReplicateClusterOut] = []
    consensus_profile: list[float] = []
    low_frequency_cluster_count = 0
    matched: dict[str, int] = {name: 0 for name in lane_names}
    missing: dict[str, list[float]] = {name: [] for name in lane_names}
    extra: dict[str, list[float]] = {name: [] for name in lane_names}
    low_frequency: dict[str, list[float]] = {name: [] for name in lane_names}
    deviation: dict[str, float] = {name: 0.0 for name in lane_names}

    for cluster_id, members in enumerate(clusters, start=1):
        lengths = [m[0] for m in members]  # 条带池有序，簇内长度天然升序
        median = _median(lengths)
        contributing = sorted({lane_names[m[1]] for m in members})
        occurrence = len(contributing)
        in_consensus = occurrence >= req.min_occurrences

        cluster_outs.append(ReplicateClusterOut(
            cluster_id=cluster_id,
            median_length=median,
            min_length=lengths[0],
            max_length=lengths[-1],
            member_count=len(members),
            occurrence_count=occurrence,
            in_consensus=in_consensus,
            lanes=contributing,
            members=[
                ReplicateClusterMemberOut(lane=lane_names[m[1]], length=m[0])
                for m in members
            ],
        ))

        if not in_consensus:
            low_frequency_cluster_count += 1
            for m in members:
                low_frequency[lane_names[m[1]]].append(m[0])
            continue

        consensus_profile.append(median)
        covered: set[str] = set()
        for name in contributing:
            # 代表条带：与中位长度最接近，并列取较短者；其余为额外条带
            own = sorted(m[0] for m in members if lane_names[m[1]] == name)
            representative = min(own, key=lambda l: (abs(l - median), l))
            own.remove(representative)
            matched[name] += 1
            deviation[name] += abs(representative - median)
            extra[name].extend(own)
            covered.add(name)
        for name in lane_names:
            if name not in covered:
                missing[name].append(median)

    # 4) 逐泳道一致性评估与稳定性评分（0~100）
    expected = len(consensus_profile)
    lane_outs: list[ReplicateLaneConsistencyOut] = []
    scores: list[float] = []
    for name in lane_names:
        hit = matched[name]
        consensus_deviation = round(deviation[name], 2)
        mean_deviation = round(consensus_deviation / hit, 2) if hit else 0.0
        missing_bands = sorted(missing[name])
        extra_bands = sorted(extra[name])
        low_frequency_bands = sorted(low_frequency[name])
        anomaly_count = (
            len(missing_bands) + len(extra_bands) + len(low_frequency_bands)
        )
        # 分母 = 共识条带数 + 额外 + 低频 = 命中数 + 异常数，恒 ≥ 1
        denominator = expected + len(extra_bands) + len(low_frequency_bands)
        score = round(100 * hit / denominator / (1 + mean_deviation), 2)
        scores.append(score)
        lane_outs.append(ReplicateLaneConsistencyOut(
            lane=name,
            band_count=len(bands_by_lane[name]),
            matched_consensus_count=hit,
            missing_count=len(missing_bands),
            missing_bands=missing_bands,
            extra_count=len(extra_bands),
            extra_bands=extra_bands,
            low_frequency_count=len(low_frequency_bands),
            low_frequency_bands=low_frequency_bands,
            anomaly_count=anomaly_count,
            consensus_deviation=consensus_deviation,
            mean_deviation=mean_deviation,
            stability_score=score,
            description=(
                f"泳道 {name}：观测条带 {len(bands_by_lane[name])} 条，"
                f"命中共识条带 {hit}/{expected}，缺失 {len(missing_bands)} 条、"
                f"额外 {len(extra_bands)} 条、低频 {len(low_frequency_bands)} 条"
                f"（异常共 {anomaly_count} 条），共识偏差 {consensus_deviation:g} bp"
                f"（平均 {mean_deviation:g} bp），稳定性评分 {score}"
            ),
        ))

    mean_score = round(sum(scores) / len(scores), 2)
    total_band_count = len(pool)
    score_reason = (
        f"全部 {total_band_count} 条观测条带按容差 {tolerance:g} bp 在 "
        f"{len(lane_names)} 条泳道间做单连接全局聚类（并列按长度与泳道序号破同分，"
        f"与输入顺序无关），共 {len(clusters)} 个簇；其中 {expected} 个簇的出现次数"
        f"（贡献泳道数）达到最少出现次数 {req.min_occurrences}，以簇中位长度构成共识谱，"
        f"其余 {low_frequency_cluster_count} 个簇为低频簇；每条泳道的异常数 = 缺失 + "
        f"额外 + 低频条带数，稳定性评分 = 100 × 命中共识条带数 ÷（共识条带数 + "
        f"额外条带数 + 低频条带数）÷（1 + 平均共识偏差），整体稳定性评分为各泳道"
        f"评分的均值 {mean_score}"
    )

    return ReplicateConsistencyResponse(
        band_match_tolerance=tolerance,
        min_occurrences=req.min_occurrences,
        lane_count=len(lane_names),
        total_band_count=total_band_count,
        cluster_count=len(clusters),
        consensus_band_count=expected,
        low_frequency_cluster_count=low_frequency_cluster_count,
        consensus_profile=consensus_profile,
        clusters=cluster_outs,
        lanes=lane_outs,
        mean_stability_score=mean_score,
        score_reason=score_reason,
    )
