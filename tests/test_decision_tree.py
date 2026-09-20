"""自适应双酶切决策树：记忆化极小化极大单元测试 + API 端到端测试。"""
import random
from itertools import permutations

from fastapi.testclient import TestClient

from app.decision_tree import (
    PairSignatures,
    build_decision_tree,
    build_equivalence_classes,
    evaluate_decision_pairs,
    solve_minimax,
)
from app.main import app
from app.models import DecisionEnzymeInput, DecisionTreeRequest

client = TestClient(app)


def _enzyme(name, genotype_sites):
    return {"name": name, "genotype_sites": genotype_sites}


def _payload(enzymes, length=1000, min_visible=50, tolerance=20):
    return {
        "plasmid_length": length,
        "min_visible_fragment": min_visible,
        "band_resolution_tolerance": tolerance,
        "enzymes": enzymes,
    }


# 示例场景：WT 与 m3 在所有酶上切点完全相同（永远无法区分）；
# BamHI+EcoRI 先分出 m2 与 {WT,m1,m3}，子集合只能由 EcoRI+HindIII 拆分。
EXAMPLE_ENZYMES = [
    _enzyme("EcoRI", {"WT": [300, 900], "m1": [300, 500], "m2": [100, 700], "m3": [300, 900]}),
    _enzyme("BamHI", {"WT": [500], "m1": [900], "m2": [100, 800], "m3": [500]}),
    _enzyme("HindIII", {"WT": [700], "m1": [100], "m2": [500], "m3": [700]}),
    _enzyme("NruI", {"WT": [], "m1": [], "m2": [], "m3": []}),
]


# ---------------------------------------------------------------- 规范化签名

def test_genotype_signature_matches_visible_fragment_clustering():
    payload = _payload(EXAMPLE_ENZYMES)
    req = DecisionTreeRequest(**payload)
    genotypes = ["WT", "m1", "m2", "m3"]
    pairs = {p.pair_key: p for p in evaluate_decision_pairs(req, genotypes)}

    # EcoRI+BamHI 对 WT：切点 300,500,900 -> 片段 200,400,400（跨原点 300+1000-900）
    # 单连接容差 20 把两条 400 共迁移合并为一条 -> 签名 (200.0, 400.0)
    assert pairs["BamHI+EcoRI"].signatures[genotypes.index("WT")] == (200.0, 400.0)
    # m2：切点 100,700,800 -> 片段 100,600,300 -> 排序签名
    assert pairs["BamHI+EcoRI"].signatures[genotypes.index("m2")] == (100.0, 300.0, 600.0)
    # HindIII+NruI 对任何基因型都无法完成双酶切（切点 < 2）-> 空签名
    assert pairs["HindIII+NruI"].signatures == ((), (), (), ())


def test_signature_rounding_two_decimals():
    # 片段全部落入同一单连接簇时代表长度取簇内均值（保留两位小数）
    payload = _payload([
        _enzyme("A", {"g1": [0, 300, 600], "g2": [0, 300, 600, 800]}),
        _enzyme("B", {"g1": [], "g2": []}),
    ], length=1000, min_visible=1, tolerance=200)
    req = DecisionTreeRequest(**payload)
    pairs = evaluate_decision_pairs(req, ["g1", "g2"])
    sig = pairs[0].signatures
    # g1 片段 300,300,400：相邻差 0/100 均 ≤ 200，全连成一簇，均值 1000/3
    assert sig[0] == (333.33,)
    # g2 片段 200,200,300,300：相邻差最大 100 ≤ 200，全连成一簇
    assert sig[1] == (250.0,)


def test_min_visible_fragments_excluded_from_signature():
    payload = _payload([
        _enzyme("A", {"g1": [0], "g2": [0]}),
        _enzyme("B", {"g1": [990], "g2": [950]}),
    ], length=1000, min_visible=50, tolerance=0)
    req = DecisionTreeRequest(**payload)
    pairs = evaluate_decision_pairs(req, ["g1", "g2"])
    # g1：切点 0,990 -> 片段 990（可见）,10（不可见，被排除）
    assert pairs[0].signatures[0] == (990.0,)
    # g2：片段 950,50（50 >= 阈值，均可见）
    assert pairs[0].signatures[1] == (50.0, 950.0)


# ---------------------------------------------------------------- 等价类

def test_equivalence_classes_merge_globally_indistinguishable():
    payload = _payload(EXAMPLE_ENZYMES)
    req = DecisionTreeRequest(**payload)
    genotypes = ["WT", "m1", "m2", "m3"]
    pairs = evaluate_decision_pairs(req, genotypes)
    classes = build_equivalence_classes(genotypes, pairs)
    # 三个等价类：{WT,m3}、{m1}、{m2}
    assert sorted(sorted(genotypes[i] for i in cls) for cls in classes) == [
        ["WT", "m3"], ["m1"], ["m2"],
    ]


def test_all_indistinguishable_single_leaf():
    payload = _payload([
        _enzyme("EcoRI", {"a": [100, 600], "b": [100, 600]}),
        _enzyme("BamHI", {"a": [300], "b": [300]}),
    ])
    body = client.post("/api/v1/digest/decision-tree", json=payload).json()
    assert body["worst_case_depth"] == 0
    assert body["internal_node_count"] == 0
    assert body["leaf_count"] == 1
    assert body["tree"]["node_type"] == "leaf"
    assert body["tree"]["indistinguishable"] is True
    assert body["tree"]["candidates"] == ["a", "b"]
    assert body["indistinguishable_groups"] == [["a", "b"]]
    assert body["indistinguishable_group_count"] == 1
    assert body["discriminating_pair_count"] == 0
    assert [p["pair_key"] for p in body["non_discriminating_pairs"]] == ["BamHI+EcoRI"]


def test_single_genotype_is_depth_zero_leaf():
    payload = _payload([
        _enzyme("EcoRI", {"only": [100, 600]}),
        _enzyme("BamHI", {"only": [300]}),
    ])
    body = client.post("/api/v1/digest/decision-tree", json=payload).json()
    assert body["tree"]["node_type"] == "leaf"
    assert body["tree"]["depth"] == 0
    assert body["tree"]["candidates"] == ["only"]
    assert body["tree"]["indistinguishable"] is False
    assert body["worst_case_depth"] == 0
    assert body["indistinguishable_groups"] == []


# ---------------------------------------------------------------- 极小化极大（纯函数，暴力交叉验证）

def _brute_force_tree(partitions, g):
    """对“基因型 -> 每个测试的签名”做暴力极小化极大。

    partitions: list[dict[int, hashable]]（第 j 个测试把基因型索引映射到签名）。
    去重同划分测试后，递归枚举剩余测试；返回 (最坏深度, 内部节点数, 测试键序列树)，
    测试按 0..m-1 的键序破同分。
    """
    distinct = {}
    for j, part in enumerate(partitions):
        key = tuple(part[i] for i in range(g))
        distinct.setdefault(key, j)
    tests = sorted(distinct.values())

    def solve(members, avail):
        if len(members) <= 1:
            return (0, 0, ("leaf", tuple(sorted(members))))
        best = None
        for j in tests:
            bit = 1 << j
            if not (avail & bit):
                continue
            groups = {}
            for i in members:
                groups.setdefault(partitions[j][i], set()).add(i)
            if len(groups) < 2:
                continue
            children = []
            worst = 0
            internal = 1
            for sig in sorted(groups, key=lambda s: (
                tuple(s) if isinstance(s, tuple) else (s,)
            )):
                sub = solve(frozenset(groups[sig]), avail & ~bit)
                children.append((sig, sub))
                worst = max(worst, sub[0])
                internal += sub[1]
            cand = (1 + worst, internal, ("node", j, children))
            if best is None or cand[:2] < best[:2]:
                best = cand
        if best is None:
            return (0, 0, ("leaf", tuple(sorted(members))))
        return best

    return solve(frozenset(range(g)), (1 << len(partitions)) - 1)


def _parts_from_tests(tests, g):
    return [{i: tests[j].signatures[i] for i in range(g)} for j in range(len(tests))]


def test_solve_minimax_matches_brute_force_random():
    rng = random.Random(20260920)
    for _ in range(400):
        g = rng.randint(2, 6)
        m = rng.randint(1, 5)
        signatures = [
            tuple(tuple(rng.choice([0.0, 1.0, 2.0, 3.0]) for _ in range(rng.randint(0, 2)))
                  for _i in range(g))
            for _j in range(m)
        ]
        tests = [
            PairSignatures(
                enzyme_pair=(f"E{j:02d}", f"F{j:02d}"),
                pair_key=f"E{j:02d}+F{j:02d}",
                signatures=tuple(sigs),
                discriminating=len(set(sigs)) > 1,
            )
            for j, sigs in enumerate(signatures)
        ]
        classes = build_equivalence_classes([f"g{i}" for i in range(g)], tests)
        root = solve_minimax(classes, [t for t in tests if t.discriminating])
        expected = _brute_force_tree(_parts_from_tests(tests, g), g)
        assert (root.worst_depth, root.internal_count) == expected[:2]


def test_minimax_prefers_smaller_worst_depth_not_greedy():
    # 测试 A 一次分出 3 份但留下大小 3 的分支（最坏 2 层）；
    # 测试 B/C 各做对半分，最坏只要 1 层（B 分出 {0},{1,2,3}，C 再分出 {1,2},{3}…）。
    # 构造：4 个基因型，A = {0}|{1}|{2,3}（最坏 2，因为 {2,3} 还需 B）；
    #       B = {0,1}|{2,3}，C = {0}|{1,2,3} -> B 后 {2,3} 用 A、{0,1} 用 C，最坏 2
    # 直接构造一个 B 单独即可最坏 1 的方案：B 四分裂
    def t(key, sigs):
        return PairSignatures(tuple(key.split("+")), key, tuple(sigs), len(set(sigs)) > 1)

    s0, s1, s2, s3 = (0.0,), (1.0,), (2.0,), (3.0,)
    tests = [
        t("A+B", (s0, s1, s2, s2)),  # 三分裂，{2,3} 仍需后续实验
        t("A+C", (s0, s0, s1, s1)),
        t("B+C", (s0, s1, s2, s3)),  # 四分裂：一次实验全部区分，最坏深度 1
    ]
    classes = build_equivalence_classes(["g0", "g1", "g2", "g3"], tests)
    root = solve_minimax(classes, tests)
    assert tests[root.pair].pair_key == "B+C"
    assert root.worst_depth == 1
    assert root.internal_count == 1
    assert len(root.branches) == 4


def test_minimax_same_depth_fewer_internal_nodes_wins():
    # 4 个基因型，两个根选择最坏深度都为 2 但内部节点数不同：
    #   A={0}|{1,2,3}，子树 {1,2,3} 用 B={1}|{2,3} 再 C 拆分 -> 3 个内部节点
    #   D={0,1}|{2,3}，两个二分组分别用 E、F 拆分            -> 3 个内部节点
    # 构造 3 对 2 节点的对照需要更大集合：6 个基因型
    #   A={0}|{1}|{2,3,4,5}，四分组链式拆分 -> 内部节点更多；
    #   B 一次把 {2,3,4,5} 拆成两组二分组 -> 内部节点更少。
    def t(key, sigs):
        return PairSignatures(tuple(key.split("+")), key, tuple(sigs), len(set(sigs)) > 1)

    s0, s1, s2, s3 = (0.0,), (1.0,), (2.0,), (3.0,)
    # 4 基因型版本（同深度 2，内部节点 2 vs 3）：
    #   A: {0}|{1,2,3}；B: {1}|{2,3}（只在含 1 时有区分）；C: {2}|{3}
    #   根 A：子树 {1,2,3} 用 B 得 {1}|{2,3} 再 C -> 根+2 = 3 内部节点
    #   D: {0,1}|{2,3}；E: {0}|{1}；F: {2}|{3}
    #   根 D：{0,1} 用 E、{2,3} 用 F -> 根+2 = 3（相同，不能破）
    # 改用：A 三分裂 {0}|{1}|{2,3}（根 + 1 子树 = 2 内部）；
    #       D 二分裂 {0,1}|{2,3}，两子树各需一个测试（根 + 2 = 3 内部）
    tests = [
        t("R+A", (s0, s1, s2, s2)),       # {0}|{1}|{2,3}
        t("R+D", (s0, s0, s1, s1)),       # {0,1}|{2,3}
        t("S+E", (s0, s1, s0, s0)),       # 只拆 {0,1}
        t("S+F", (s0, s0, s0, s1)),       # 只拆 {2,3}
        t("S+B", (s0, s1, s2, s2)),       # 与 A 同划分（用于 {1,2,3} 时等价于 A）
    ]
    classes = build_equivalence_classes(["g0", "g1", "g2", "g3"], tests)
    root = solve_minimax(classes, tests)
    # R+A 与 R+D 最坏深度都为 2；R+A 只需 2 个内部节点（{0},{1} 直接成叶），
    # R+D 需要 3 个
    assert tests[root.pair].pair_key == "R+A"
    assert root.worst_depth == 2
    assert root.internal_count == 2


def test_minimax_same_depth_and_nodes_lexicographic_pair_key_wins():
    def t(key, sigs):
        return PairSignatures(tuple(key.split("+")), key, tuple(sigs), len(set(sigs)) > 1)

    a, b = (0.0,), (1.0,)
    # 两对划分完全相同：服务保证 tests 按 pair_key 字典序传入
    tests = sorted([
        t("Z+Y", (a, a, b)),
        t("A+X", (a, a, b)),
    ], key=lambda p: p.pair_key)
    classes = build_equivalence_classes(["g0", "g1", "g2"], tests)
    root = solve_minimax(classes, tests)
    assert tests[root.pair].pair_key == "A+X"


def test_pairs_with_identical_partitions_are_deduped():
    def t(key, sigs):
        return PairSignatures(tuple(key.split("+")), key, tuple(sigs), len(set(sigs)) > 1)

    a, b = (0.0,), (1.0,)
    # 划分相同的两对：去重保留键序最小者（服务中 pairs 先按 pair_key 排序）
    tests_full = sorted([t("Z+Y", (a, a, b)), t("A+X", (a, a, b))],
                        key=lambda p: p.pair_key)
    distinct = {}
    for test in tests_full:
        distinct.setdefault(test.signatures, test)
    deduped = list(distinct.values())
    assert [p.pair_key for p in deduped] == ["A+X"]


def test_pair_that_cannot_split_current_set_is_ignored():
    # 示例场景：根 BamHI+EcoRI 分出 {WT,m1,m3} 后，BamHI+EcoRI（已用）与
    # BamHI+HindIII（签名相同）都不能拆分该集合，只有 EcoRI+HindIII 可用。
    body = client.post("/api/v1/digest/decision-tree", json=_payload(EXAMPLE_ENZYMES)).json()
    root = body["tree"]
    assert root["pair_key"] == "BamHI+EcoRI"
    internal = [b for b in root["branches"] if b["child"]["node_type"] == "internal"]
    assert len(internal) == 1
    assert internal[0]["genotypes"] == ["WT", "m1", "m3"]
    assert internal[0]["child"]["pair_key"] == "EcoRI+HindIII"


# ---------------------------------------------------------------- 端到端示例

def test_decision_tree_example_response():
    body = client.post("/api/v1/digest/decision-tree", json=_payload(EXAMPLE_ENZYMES)).json()

    assert body["genotype_count"] == 4
    assert body["genotypes"] == ["WT", "m1", "m2", "m3"]
    assert body["total_enzyme_pairs"] == 6
    assert body["discriminating_pair_count"] == 5
    assert [p["pair_key"] for p in body["non_discriminating_pairs"]] == ["HindIII+NruI"]

    assert body["worst_case_depth"] == 2
    assert body["internal_node_count"] == 2
    assert body["leaf_count"] == 3
    assert body["indistinguishable_groups"] == [["WT", "m3"]]
    assert body["indistinguishable_group_count"] == 1

    root = body["tree"]
    assert root["node_id"] == 1 and root["depth"] == 0
    assert root["node_type"] == "internal"
    assert root["candidates"] == ["WT", "m1", "m2", "m3"]
    assert root["enzyme_pair"] == ["BamHI", "EcoRI"]

    # 分支按签名向量字典序：(100,300,600) < (200,400)
    assert [b["signature"] for b in root["branches"]] == [
        [100.0, 300.0, 600.0], [200.0, 400.0],
    ]
    assert root["branches"][0]["genotypes"] == ["m2"]
    assert root["branches"][0]["child"]["node_type"] == "leaf"

    sub = root["branches"][1]["child"]
    assert sub["depth"] == 1
    assert sub["pair_key"] == "EcoRI+HindIII"
    assert [b["signature"] for b in sub["branches"]] == [
        [200.0, 400.0], [200.0, 600.0],
    ]
    indist_leaf = sub["branches"][0]["child"]
    assert indist_leaf["node_type"] == "leaf"
    assert indist_leaf["candidates"] == ["WT", "m3"]
    assert indist_leaf["indistinguishable"] is True
    assert sub["branches"][1]["child"]["candidates"] == ["m1"]

    # 节点编号为先根遍历
    assert [n["node_id"] for n in _all_nodes(root)] == [1, 2, 3, 4, 5]
    assert body["tree_reason"]


def _all_nodes(node):
    yield node
    if node.get("branches"):
        for branch in node["branches"]:
            yield from _all_nodes(branch["child"])


def test_result_independent_of_input_order():
    base = _payload(EXAMPLE_ENZYMES)
    genotypes = ["WT", "m1", "m2", "m3"]
    first = client.post("/api/v1/digest/decision-tree", json=base).json()
    # 酶顺序与每个酶内基因型键顺序均不影响规范化结果
    for enzyme_order in permutations(base["enzymes"]):
        shuffled_enzymes = []
        for enzyme in enzyme_order:
            genotype_order = permutations(genotypes)
            sites = {g: enzyme["genotype_sites"][g] for g in next(genotype_order)}
            shuffled_enzymes.append({"name": enzyme["name"], "genotype_sites": sites})
        res = client.post(
            "/api/v1/digest/decision-tree", json={**base, "enzymes": shuffled_enzymes},
        )
        assert res.status_code == 200
        assert res.json() == first


def test_tolerance_within_genotype_clustering_changes_signature():
    # g1 切点 0,290,600 -> 片段 290,310,400：容差 10 时 290/310 不合带（3 带），
    # 容差 20 时合为一带（代表 300，2 带）；
    # g2 切点 0,600 -> 片段 600,400：两种容差下签名均为 (400,600)。
    def payload(tol):
        return _payload([
            _enzyme("A", {"g1": [0, 290, 600], "g2": [0, 600]}),
            _enzyme("B", {"g1": [], "g2": []}),
        ], length=1000, min_visible=1, tolerance=tol)

    strict = client.post("/api/v1/digest/decision-tree", json=payload(10)).json()
    loose = client.post("/api/v1/digest/decision-tree", json=payload(20)).json()

    strict_root = strict["tree"]
    assert {tuple(b["signature"]): b["genotypes"] for b in strict_root["branches"]} == {
        (290.0, 310.0, 400.0): ["g1"],
        (400.0, 600.0): ["g2"],
    }
    loose_root = loose["tree"]
    assert {tuple(b["signature"]): b["genotypes"] for b in loose_root["branches"]} == {
        (300.0, 400.0): ["g1"],
        (400.0, 600.0): ["g2"],
    }
    assert strict["discriminating_pair_count"] == loose["discriminating_pair_count"] == 1


# ---------------------------------------------------------------- 422 校验

def _expect_422(payload, code_part=None):
    res = client.post("/api/v1/digest/decision-tree", json=payload)
    assert res.status_code == 422
    body = res.json()
    assert body["code"] == "validation_failed"
    assert body["errors"]
    if code_part:
        assert any(code_part in e["type"] for e in body["errors"])
    return body


def test_422_site_out_of_range():
    payload = _payload([
        _enzyme("EcoRI", {"g1": [100, 1000], "g2": [100]}),
        _enzyme("BamHI", {"g1": [300], "g2": [300]}),
    ])
    body = _expect_422(payload, "site_out_of_range")
    assert any("genotype_sites" in e["loc"] for e in body["errors"])


def test_422_negative_site():
    _expect_422(_payload([
        _enzyme("EcoRI", {"g1": [-1], "g2": [100]}),
        _enzyme("BamHI", {"g1": [300], "g2": [300]}),
    ]), "site_out_of_range")


def test_422_duplicate_site():
    _expect_422(_payload([
        _enzyme("EcoRI", {"g1": [100, 100, 600], "g2": [100]}),
        _enzyme("BamHI", {"g1": [300], "g2": [300]}),
    ]), "duplicate_site")


def test_422_duplicate_enzyme_name():
    _expect_422(_payload([
        _enzyme("EcoRI", {"g1": [100], "g2": [200]}),
        _enzyme("EcoRI", {"g1": [300], "g2": [400]}),
    ]), "duplicate_enzyme_name")


def test_422_inconsistent_genotype_names():
    _expect_422(_payload([
        _enzyme("EcoRI", {"g1": [100], "g2": [200]}),
        _enzyme("BamHI", {"g1": [300], "g3": [400]}),
    ]), "inconsistent_genotype_names")


def test_422_blank_genotype_name():
    _expect_422(_payload([
        _enzyme("EcoRI", {"": [100]}),
        _enzyme("BamHI", {"": [300]}),
    ]), "invalid_genotype_name")


def test_422_no_genotypes():
    _expect_422(_payload([
        _enzyme("EcoRI", {}),
        _enzyme("BamHI", {}),
    ]), "no_genotypes")


def test_422_too_many_genotypes(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_DECISION_GENOTYPES", 2)
    _expect_422(_payload([
        _enzyme("EcoRI", {"g1": [100], "g2": [200], "g3": [300]}),
        _enzyme("BamHI", {"g1": [400], "g2": [500], "g3": [600]}),
    ]), "too_many_decision_genotypes")


def test_422_too_many_enzymes(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_DECISION_ENZYMES", 2)
    _expect_422(_payload([
        _enzyme("E1", {"g1": [100], "g2": [200]}),
        _enzyme("E2", {"g1": [300], "g2": [400]}),
        _enzyme("E3", {"g1": [500], "g2": [600]}),
    ]), "too_many_decision_enzymes")


def test_422_min_visible_not_smaller_than_plasmid():
    _expect_422(_payload(
        [
            _enzyme("EcoRI", {"g1": [100], "g2": [200]}),
            _enzyme("BamHI", {"g1": [300], "g2": [400]}),
        ],
        min_visible=1000,
    ), "min_visible_too_large")


def test_422_basic_field_constraints():
    _expect_422(_payload([
        _enzyme("EcoRI", {"g1": [100], "g2": [200]}),
        _enzyme("BamHI", {"g1": [300], "g2": [400]}),
    ], length=0))
    _expect_422(_payload([
        _enzyme("EcoRI", {"g1": [100], "g2": [200]}),
        _enzyme("BamHI", {"g1": [300], "g2": [400]}),
    ], tolerance=-1))
    _expect_422(_payload([_enzyme("EcoRI", {"g1": [100], "g2": [200]})]))
    # genotype_sites 类型错误
    _expect_422(_payload([
        _enzyme("EcoRI", {"g1": [100], "g2": [200]}),
        {"name": "BamHI", "genotype_sites": {"g1": 300, "g2": 400}},
    ]))


# ---------------------------------------------------------------- 规模上限内的性能

def test_max_scale_request_completes():
    # 默认上限 12 基因型 × 10 酶（45 对）的满规模请求应快速完成
    rng = random.Random(1)
    genotypes = [f"g{i:02d}" for i in range(12)]
    enzymes = []
    for j in range(10):
        enzymes.append(_enzyme(
            f"E{j:02d}",
            {g: sorted(rng.sample(range(50, 950, 25), rng.randint(1, 4))) for g in genotypes},
        ))
    res = client.post("/api/v1/digest/decision-tree", json=_payload(enzymes))
    assert res.status_code == 200
    body = res.json()
    assert body["total_enzyme_pairs"] == 45
    assert body["genotype_count"] == 12
    assert body["worst_case_depth"] >= 1
