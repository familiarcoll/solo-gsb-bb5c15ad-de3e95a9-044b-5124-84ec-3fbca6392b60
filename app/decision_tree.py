"""核心算法：自适应双酶切决策树（记忆化极小化极大，精确最小化最坏实验次数）。

所有函数均为纯函数，不依赖任何外部状态；结果只取决于切点数据本身，与酶、
基因型的提交顺序无关。

规范化条带签名
--------------
对每一对不同酶、每一个基因型，复用 :mod:`app.digest` 的圆环双酶切规则合并
两酶切点并按圆环相邻切点计算片段（跨原点片段长度为 p0 + L − p_last），保留
长度 ≥ 最小可见片段的可见片段，再复用 :mod:`app.mixture` 的单连接聚带规则
（相邻片段长度差 ≤ 分辨容差归入同一簇，代表长度取簇内均值并保留两位小数），
得到该基因型在该酶对上**规范化的条带签名**（代表长度元组，升序）。某基因型
在某酶对上有效切点不足 2 个时没有任何片段，签名为空元组。

不可区分组合并
--------------
一个基因型在**全部酶对**上的签名构成其签名向量；两个基因型的签名向量完全
相同，当且仅当没有任何酶对能区分它们（任一酶对上二者都落入同一签名），即
永远无法区分。决策树直接在这些等价类上构建。

记忆化极小化极大
----------------
状态为 (候选等价类集合, 剩余酶对集合)。内部节点选择一个剩余酶对做实验，按
观测到的签名把候选类拆成若干分支，递归求解；节点代价为
``1 + max(各分支最坏深度)``。搜索依次：

  1. **最坏深度最小**（极小化极大，精确最优，非贪心）；
  2. 最坏深度相同时，**内部节点总数更少**（多分支拆分优先）；
  3. 仍相同时，**酶对 pair_key 字典序更小**（酶对按键序枚举，首个最优保留）。

不能拆分当前候选集合的酶对直接忽略；同一基因型划分完全相同的酶对在搜索中
等价，只保留 pair_key 最小的代表。状态做记忆化，结果唯一且与输入顺序无关。
"""
from dataclasses import dataclass, field
from itertools import combinations

from app.digest import circular_fragments, merge_cut_points, pair_key
from app.mixture import cluster_predicted_bands
from app.models import (
    DecisionBranchOut,
    DecisionPairOut,
    DecisionTreeNodeOut,
    DecisionTreeRequest,
    DecisionTreeResponse,
)


# --------------------------------------------------------------- 规范化签名

def genotype_signature(site_map: dict[str, dict[str, list[int]]],
                       enzyme_a: str, enzyme_b: str, genotype: str,
                       plasmid_length: int, min_visible: int,
                       tolerance: int) -> tuple[float, ...]:
    """单个基因型在单个酶对上的规范化条带签名。"""
    cuts = merge_cut_points(site_map, enzyme_a, enzyme_b, genotype)
    if len(cuts) < 2:
        return ()
    fragments = circular_fragments(cuts, plasmid_length)
    lengths = [f.length for f in fragments if f.length >= min_visible]
    return tuple(cluster_predicted_bands(lengths, tolerance))


@dataclass(frozen=True)
class PairSignatures:
    """一个酶对在全部基因型（按基因型字典序）上的规范化签名。"""

    enzyme_pair: tuple[str, str]
    pair_key: str
    signatures: tuple[tuple[float, ...], ...]
    discriminating: bool = field(compare=False)


def evaluate_decision_pairs(req: DecisionTreeRequest,
                            genotypes: list[str]) -> list[PairSignatures]:
    """枚举全部酶对并计算逐基因型签名；酶对按 pair_key 字典序排列。"""
    site_map = {
        enzyme.name: {g: sorted(enzyme.genotype_sites[g]) for g in genotypes}
        for enzyme in req.enzymes
    }
    names = sorted(enzyme.name for enzyme in req.enzymes)

    pairs: list[PairSignatures] = []
    for a, b in combinations(names, 2):
        _, key = pair_key(a, b)
        signatures = tuple(
            genotype_signature(
                site_map, a, b, g,
                req.plasmid_length,
                req.min_visible_fragment,
                req.band_resolution_tolerance,
            )
            for g in genotypes
        )
        pairs.append(PairSignatures(
            enzyme_pair=(a, b),
            pair_key=key,
            signatures=signatures,
            discriminating=len(set(signatures)) > 1,
        ))
    return pairs


# --------------------------------------------------------------- 等价类与搜索结构

@dataclass
class _Node:
    """极小化极大搜索的内部树节点（成员为等价类索引）。"""

    members: tuple[int, ...]
    pair: int | None = None                # 选中的酶对（tests 列表索引）；叶节点为 None
    branches: list[tuple[tuple[float, ...], "_Node"]] | None = None
    is_leaf: bool = False
    worst_depth: int = 0                   # 以该节点为根的最坏深度（叶为 0）
    internal_count: int = 0                # 以该节点为根的内部节点总数（叶为 0）


def build_equivalence_classes(
        genotypes: list[str],
        pairs: list[PairSignatures]) -> list[list[int]]:
    """按全部酶对上的签名向量合并永远无法区分的基因型，返回等价类（基因型索引）。

    各等价类按其基因型名元组字典序排列，类内基因型索引亦按名字典序排列。
    """
    vectors = [
        tuple(p.signatures[i] for p in pairs) for i in range(len(genotypes))
    ]
    groups: dict[tuple, list[int]] = {}
    for i, vector in enumerate(vectors):
        groups.setdefault(vector, []).append(i)
    classes = [
        sorted(indices, key=lambda i: genotypes[i]) for indices in groups.values()
    ]
    classes.sort(key=lambda cls: tuple(genotypes[i] for i in cls))
    return classes


def solve_minimax(classes: list[list[int]],
                  tests: list[PairSignatures]) -> _Node:
    """在等价类上做记忆化极小化极大搜索，返回最优（子）树根节点。

    tests 为去重后的可拆分酶对（按 pair_key 字典序）；状态为
    (候选等价类位掩码, 剩余酶对位掩码)。
    """
    memo: dict[tuple[int, int], _Node] = {}
    full_pairs = (1 << len(tests)) - 1
    # 等价类内签名一致，取类内首个基因型的签名作为该类在各酶对上的签名
    class_signatures = [
        [test.signatures[cls[0]] for cls in classes] for test in tests
    ]

    def solve(class_mask: int, pair_mask: int) -> _Node:
        key = (class_mask, pair_mask)
        cached = memo.get(key)
        if cached is not None:
            return cached

        members = tuple(
            i for i in range(len(classes)) if (class_mask >> i) & 1
        )
        if len(members) <= 1:
            leaf = _Node(members=members, is_leaf=True)
            memo[key] = leaf
            return leaf

        best: tuple[int, int, _Node] | None = None  # (最坏深度, 内部节点数, 节点)
        for j, test in enumerate(tests):
            bit = 1 << j
            if not (pair_mask & bit):
                continue
            grouped: dict[tuple[float, ...], int] = {}
            for i in members:
                sig = class_signatures[j][i]
                grouped[sig] = grouped.get(sig, 0) | (1 << i)
            if len(grouped) < 2:
                # 该酶对不能拆分当前集合：忽略
                continue

            child_avail = pair_mask & ~bit
            children: list[tuple[tuple[float, ...], _Node]] = []
            worst = 0
            internal = 1
            for sig in sorted(grouped):
                child = solve(grouped[sig], child_avail)
                children.append((sig, child))
                worst = max(worst, child.worst_depth)
                internal += child.internal_count
            node = _Node(
                members=members, pair=j, branches=children,
                worst_depth=1 + worst, internal_count=internal,
            )
            candidate = (node.worst_depth, node.internal_count, node)
            if best is None or candidate[:2] < best[:2]:
                # tests 按 pair_key 字典序枚举：同深度同内部节点数时保留最先遇到者
                best = candidate
            if best[0] == 1:
                # 一次实验即把全部候选类分入单类分支：理论最优，键序最小者即最终解
                break

        if best is None:
            # 理论上不可达：不同等价类必有至少一个可拆分酶对；防御性返回叶节点
            leaf = _Node(members=members, is_leaf=True)
            memo[key] = leaf
            return leaf

        memo[key] = best[2]
        return best[2]

    initial_mask = (1 << len(classes)) - 1
    return solve(initial_mask, full_pairs)


# --------------------------------------------------------------- 响应组装

def _serialize(node: _Node, depth: int, counter: list[int],
               classes: list[list[int]], genotypes: list[str],
               tests: list[PairSignatures]) -> DecisionTreeNodeOut:
    counter[0] += 1
    node_id = counter[0]
    candidates = sorted(
        (genotypes[i] for cls in node.members for i in classes[cls]),
    )

    if node.is_leaf:
        return DecisionTreeNodeOut(
            node_id=node_id,
            depth=depth,
            node_type="leaf",
            candidates=candidates,
            indistinguishable=len(candidates) > 1,
        )

    test = tests[node.pair]  # type: ignore[arg-type]
    branches = [
        DecisionBranchOut(
            signature=list(sig),
            genotypes=sorted(
                genotypes[i] for cls in child.members for i in classes[cls]
            ),
            child=_serialize(child, depth + 1, counter, classes, genotypes, tests),
        )
        for sig, child in node.branches  # type: ignore[union-attr]
    ]
    return DecisionTreeNodeOut(
        node_id=node_id,
        depth=depth,
        node_type="internal",
        candidates=candidates,
        pair_key=test.pair_key,
        enzyme_pair=list(test.enzyme_pair),
        branches=branches,
    )


def _count_leaves(node: DecisionTreeNodeOut) -> int:
    if node.node_type == "leaf":
        return 1
    return sum(_count_leaves(b.child) for b in node.branches)  # type: ignore[union-attr]


def build_decision_tree(req: DecisionTreeRequest) -> DecisionTreeResponse:
    genotypes = sorted({g for e in req.enzymes for g in e.genotype_sites})
    pairs = evaluate_decision_pairs(req, genotypes)

    non_discriminating = [p for p in pairs if not p.discriminating]
    discriminating = [p for p in pairs if p.discriminating]

    # 划分完全相同的酶对在搜索中等价：只保留 pair_key 字典序最小的代表
    distinct: dict[tuple[tuple[float, ...], ...], PairSignatures] = {}
    for test in discriminating:  # pairs 已按 pair_key 排序，首次出现即键序最小
        distinct.setdefault(test.signatures, test)
    tests = list(distinct.values())

    classes = build_equivalence_classes(genotypes, pairs)
    root = solve_minimax(classes, tests)

    counter = [0]
    tree = _serialize(root, 0, counter, classes, genotypes, tests)
    worst_depth = root.worst_depth
    internal_count = root.internal_count
    leaf_count = _count_leaves(tree)

    indistinguishable_groups = [
        sorted(genotypes[i] for i in cls)
        for cls in classes if len(cls) > 1
    ]

    reason = (
        f"共枚举 {len(pairs)} 个酶对：{len(discriminating)} 个酶对至少能拆分一对基因型，"
        f"{len(non_discriminating)} 个酶对对全部基因型签名相同（忽略）；"
        f"{len(discriminating) - len(tests)} 个酶对的基因型划分与更小 pair_key 的酶对完全"
        "等价（搜索去重）。每个基因型在每个酶对上按圆环双酶切片段经分辨容差单连接聚带、"
        "簇内均值保留两位小数生成规范化条带签名；在全部酶对上签名向量相同的基因型合并为"
        f"{len(indistinguishable_groups)} 个不可区分组。随后以“候选等价类集合 + 剩余酶对”"
        "为状态做记忆化极小化极大搜索，精确最小化根到叶的最坏深度；同深度依次取内部节点"
        "更少、酶对 pair_key 字典序更小的方案，且不选择不能拆分当前集合的酶对，"
        f"得最坏实验次数 {worst_depth}、内部节点 {internal_count} 个。全部规范化只依赖切点"
        "数据，结果与酶和基因型的输入顺序无关。"
    )

    return DecisionTreeResponse(
        plasmid_length=req.plasmid_length,
        min_visible_fragment=req.min_visible_fragment,
        band_resolution_tolerance=req.band_resolution_tolerance,
        genotype_count=len(genotypes),
        genotypes=genotypes,
        total_enzyme_pairs=len(pairs),
        discriminating_pair_count=len(discriminating),
        non_discriminating_pairs=[
            DecisionPairOut(enzyme_pair=list(p.enzyme_pair), pair_key=p.pair_key)
            for p in non_discriminating
        ],
        tree=tree,
        worst_case_depth=worst_depth,
        internal_node_count=internal_count,
        leaf_count=leaf_count,
        indistinguishable_group_count=len(indistinguishable_groups),
        indistinguishable_groups=indistinguishable_groups,
        tree_reason=reason,
    )
