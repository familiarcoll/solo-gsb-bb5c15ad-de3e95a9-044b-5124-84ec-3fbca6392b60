# 环状 DNA 双酶切鉴别服务

基于 FastAPI 的圆形（环状）DNA 双酶切鉴别服务，提供八个互补接口：

- **方案推荐**：客户端提交质粒长度、参考型与变异型的候选酶切位点、最小可见片段长度和条带分辨容差，服务枚举所有不同酶对，按圆环相邻切点计算两种样本的片段长度，将容差内的片段聚为同一条电泳条带，并推荐能够区分两种样本的酶对。
- **样本鉴别**：客户端提交同一质粒参考型、变异型的预测条带长度和一次凝胶观测长度列表，服务在给定容差内做确定性的一对一最小代价匹配，分别计算两种样本的缺失条带、额外条带与总误差，判定观测样本为参考型、变异型或无法区分。
- **最小酶对面板**：客户端提交质粒长度、可见阈值、分辨容差及酶列表（每种酶给出参考型切点和按变异名映射的变异型切点），服务复用同一套圆环酶切与聚带规则评估每个酶对对每个变异型的可区分条带，并**精确**求出覆盖全部变异型的最小酶对面板（非贪心近似）。
- **混合样本解析**：客户端提交若干**命名候选成分**各自的预测条带、一次混合样本观测条带、匹配容差与成分数上限，服务枚举不超过上限的全部候选子集，先把每个子集的预测条带按容差做单连接聚带（共迁移条带只出现一次、代表长度取均值），再复用鉴别接口的一对一最小代价匹配，按未匹配条带数、总误差、成分数、候选名列表字典序唯一选定最吻合的成分组合。
- **批量盲样分配**：客户端提交若干**命名预测谱型**与若干**命名观测泳道**、条带匹配容差和单条分配允许的最大未匹配条带数，服务先复用鉴别接口的一对一条带匹配评估每条“谱型-泳道”边（未匹配数超阈值的边不可用），再在可用边上做**全局一对一分配**，依次最大化已分配泳道数、最小化未匹配总数、最小化总误差，并以泳道名排序后的谱型名序列字典序破同分，返回逐项匹配依据、未分配泳道与未使用谱型。
- **自适应决策树**：客户端提交质粒长度、可见阈值、分辨容差及各基因型在各酶上的切点，服务为每个酶对 × 每个基因型复用同一套圆环酶切与聚带规则生成**规范化条带签名**，先按全部酶对上的签名向量合并永远无法区分的基因型，再以“候选等价类集合 + 剩余酶对”为状态做**记忆化极小化极大搜索**，精确生成最小化最坏实验次数的自适应决策树（同深度依次取内部节点更少、酶对键字典序更小的方案，忽略不能拆分当前集合的酶对）。
- **部分酶切失败反推**：客户端提交质粒长度、按圆环顺序排列的理论片段（每段长度及其后边界名称）、一次异常凝胶的观测条带、匹配容差与最大失败边界数，服务枚举失败边界数不超过上限的**全部边界子集**；失败（未切开）边界会合并其相邻片段（正确处理连续失败与首尾**跨原点**合并），再复用鉴别接口的一对一条带匹配，按未匹配条带数、总误差、失败边界数、边界名序列依次选出**唯一诊断**，返回合并后预测条带、逐条匹配与失败边界。
- **凝胶重复实验一致性评估**：客户端提交同一样本的若干**命名观测泳道**、条带匹配容差与最少出现次数，服务把所有泳道的条带汇入全局条带池做**确定性的单连接聚类**（并列按长度与泳道序号破同分，结果与输入顺序无关），以出现次数（贡献泳道数）达到阈值的簇的**中位长度**构成共识谱，标记每条泳道的**缺失、额外及低频条带**，并按共识偏差与异常数给出 0~100 的**稳定性评分**。

**核心功能为纯本地数学计算，不依赖数据库、缓存或任何外部服务。**

---

## 单条命令启动

需要 Docker（含 Compose 插件，Docker 20.10+ / Docker Desktop 自带）：

```bash
docker compose up --build -d
```

服务启动后监听 **`http://localhost:8000`**。

停止：

```bash
docker compose down
```

### 不使用 Docker（本地 Python 3.11+）

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

---

## API 地址与端口

| 用途 | 方法 | 地址 |
| --- | --- | --- |
| 健康检查 | GET | `http://localhost:8000/health` |
| 双酶切方案推荐 | POST | `http://localhost:8000/api/v1/digest/recommend` |
| 凝胶条带样本鉴别 | POST | `http://localhost:8000/api/v1/digest/identify` |
| 最小酶对面板 | POST | `http://localhost:8000/api/v1/digest/panel` |
| 混合样本解析 | POST | `http://localhost:8000/api/v1/digest/mixture` |
| 批量盲样分配 | POST | `http://localhost:8000/api/v1/digest/batch-identify` |
| 自适应决策树 | POST | `http://localhost:8000/api/v1/digest/decision-tree` |
| 部分酶切失败反推 | POST | `http://localhost:8000/api/v1/digest/partial-failure` |
| 凝胶重复实验一致性评估 | POST | `http://localhost:8000/api/v1/digest/replicate-consistency` |
| 交互式 API 文档（Swagger） | GET | `http://localhost:8000/docs` |
| OpenAPI Schema | GET | `http://localhost:8000/openapi.json` |

健康检查示例：

```bash
curl http://localhost:8000/health
# {"status":"ok","service":"circular-dna-digest","version":"1.7.0"}
```

---

## 必要配置

所有配置均有默认值，**开箱即用，无需任何外部服务**。可通过环境变量覆盖：

| 环境变量 | 默认值 | 说明 | Compose 中如何改 |
| --- | --- | --- | --- |
| `APP_HOST` | `0.0.0.0` | 服务监听地址 | 编辑 `docker-compose.yml` 的 `environment` |
| `APP_PORT` | `8000` | 容器内监听端口 | 编辑 `docker-compose.yml` 的 `environment` |
| `HOST_PORT` | `8000` | 宿主机映射端口 | `HOST_PORT=9000 docker compose up --build -d` |
| `MAX_ENZYMES` | `64` | 单次请求允许的最大酶数量（酶对数为 C(n,2)） | 编辑 `environment` |
| `MAX_BANDS` | `512` | identify / mixture / batch-identify 接口每个条带列表允许的最大条带数 | 编辑 `environment` |
| `MAX_VARIANTS` | `16` | panel 接口允许的最大变异型数量（最小面板为精确集合覆盖，规模随变异数指数增长） | 编辑 `environment` |
| `MAX_MIXTURE_CANDIDATES` | `12` | mixture 接口允许的最大候选成分数（需枚举至多该规模的全部子集） | 编辑 `environment` |
| `MAX_BATCH_PROFILES` | `32` | batch-identify 接口允许的最大预测谱型数（候选边数为谱型数 × 泳道数） | 编辑 `environment` |
| `MAX_BATCH_LANES` | `32` | batch-identify 接口允许的最大观测泳道数（全局一对一分配为最小费用最大流） | 编辑 `environment` |
| `MAX_DECISION_GENOTYPES` | `12` | decision-tree 接口允许的最大基因型数（极小化极大搜索状态随基因型数增长） | 编辑 `environment` |
| `MAX_DECISION_ENZYMES` | `10` | decision-tree 接口允许的最大酶数量（酶对数为 C(n,2)，独立于 `MAX_ENZYMES`） | 编辑 `environment` |
| `MAX_PARTIAL_FAILURE_BOUNDARIES` | `16` | partial-failure 接口允许的最大理论片段（边界）数（需枚举至多 `max_failed_boundaries` 规模的全部边界子集，最坏 2ⁿ） | 编辑 `environment` |
| `MAX_REPLICATE_LANES` | `32` | replicate-consistency 接口允许的最大观测泳道数（全局聚类规模为泳道数 × 每泳道条带数） | 编辑 `environment` |

---

## 请求与响应

### 双酶切方案推荐（POST `/api/v1/digest/recommend`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `plasmid_length` | int | > 0，环状质粒总长度（bp），切点合法范围为 `0 ~ plasmid_length-1` |
| `min_visible_fragment` | int | ≥ 0 且 < 质粒长度；短于此长度的片段不参与条带聚类 |
| `band_resolution_tolerance` | int | ≥ 0，条带分辨容差（bp） |
| `enzymes` | array | 至少 2 种酶，酶名称不可重复 |
| `enzymes[].name` | string | 酶名称，如 `EcoRI` |
| `enzymes[].reference_sites` | int[] | 参考型候选切点（0-based，bp） |
| `enzymes[].variant_sites` | int[] | 变异型候选切点（0-based，bp） |

请求示例见 [`examples/request.json`](examples/request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/recommend \
  -H 'Content-Type: application/json' \
  -d @examples/request.json
```

#### 算法说明（推荐）

1. **枚举酶对**：对所有不同酶的两两组合（酶名按字典序规范化）。
2. **合并切点**：每种样本合并两酶切点；两种酶切在同一位置时合并为一个切点并记录来源酶（不视为重复位点）。
3. **圆环片段**：切点升序排列为 p0 < p1 < … < p_{n-1}，相邻切点间形成片段；末段 `(p_{n-1}, p0)` 跨过坐标原点，长度为 `p0 + L - p_{n-1}`，响应中以 `crosses_origin: true` 标记。
4. **可见性过滤**：长度 < `min_visible_fragment` 的片段不参与成带（计数为 `hidden_fragment_count`）。
5. **条带聚类**：两样本可见片段合并并按长度排序，相邻片段长度差 ≤ `band_resolution_tolerance` 归入同一条带（单连接聚类），超过则断开成新带。
6. **可区分条带**：条带内参考型与变异型的片段数量不同（`reference_only` / `variant_only` / `count_mismatch`）。
7. **推荐排序**：可区分条带数降序 → 总条带数升序 → 酶对名称字典序升序。

#### 响应要点（推荐）

- `recommendations`：推荐酶对，包含
  - `samples.reference` / `samples.variant`：合并后的**切点**（含来源酶）、全部**片段**及跨原点标记、可见/不可见片段计数；
  - `bands`：聚类后的**条带**（代表长度、长度范围、两样本各自片段数、成员片段、是否可区分）；
  - `differences`：**差异依据**（中文说明 + 机器可读类型）。
- `rejected_pairs`：**任一样本有效切点不足 2 个**的酶对（含原因与各自切点数量）。
- `non_discriminating_pairs`：切点足够但两样本条带完全一致、无法鉴别的酶对。
- `total_enzyme_pairs` / `recommended_count` / `rejected_pair_count` / `non_discriminating_pair_count`：汇总计数。

### 凝胶条带样本鉴别（POST `/api/v1/digest/identify`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `reference_bands` | number[] | 参考型预测条带长度（bp），每条 ≥ 0，可为空列表 |
| `variant_bands` | number[] | 变异型预测条带长度（bp），每条 ≥ 0，可为空列表 |
| `observed_bands` | number[] | 一次凝胶观测到的条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `band_match_tolerance` | number | ≥ 0，条带匹配容差（bp） |

请求示例见 [`examples/identify-request.json`](examples/identify-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/identify \
  -H 'Content-Type: application/json' \
  -d @examples/identify-request.json
```

#### 匹配与判定规则

1. **一对一最小代价匹配**：对每种样本，将其预测条带与观测条带分别按长度升序排列后做一对一配对；仅当 `|预测长度 − 观测长度| ≤ band_match_tolerance` 时允许配对。目标为先**最大化匹配对数**，再**最小化总误差**（各配对长度差之和）。一维绝对差代价满足 Monge 性质，最优解必可取不交叉匹配，服务用序列比对式动态规划在 O(n·m) 内求得；并列时按固定优先级选择（等效于较短条带优先配对），结果完全确定且与列表提交顺序无关。
2. **缺失 / 额外条带**：未匹配到观测条带的预测条带为**缺失条带**（`missing_bands`），未匹配到任何预测条带的观测条带为**额外条带**（`extra_bands`）。
3. **判定**：比较两种样本的未匹配条带数（缺失 + 额外），少者判为对应样本（`reference` / `variant`）；并列时总误差小者胜出；两者完全一致时判定为 `ambiguous`（无法区分）。

#### 响应要点（鉴别）

- `samples.reference` / `samples.variant`：每种样本的
  - `matched_count` / `missing_count` / `extra_count` / `total_error`：匹配对数、缺失数、额外数与总误差（bp，保留两位小数）；
  - `matches`：**逐条匹配依据**（预测长度、观测长度、长度差及中文说明），按预测长度升序；
  - `missing_bands` / `extra_bands`：缺失、额外条带长度列表（升序）。
- `verdict`：`reference` / `variant` / `ambiguous`；`verdict_reason` 为中文判定依据，可完全由响应中的计数与总误差复算。

### 最小酶对面板（POST `/api/v1/digest/panel`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `plasmid_length` | int | > 0，环状质粒总长度（bp），切点合法范围为 `0 ~ plasmid_length-1` |
| `min_visible_fragment` | int | ≥ 0 且 < 质粒长度；短于此长度的片段不参与条带聚类 |
| `band_resolution_tolerance` | int | ≥ 0，条带分辨容差（bp） |
| `enzymes` | array | 至少 2 种酶，酶名称不可重复；所有酶的变异名集合必须完全一致 |
| `enzymes[].name` | string | 酶名称，如 `EcoRI` |
| `enzymes[].reference_sites` | int[] | 参考型候选切点（0-based，bp） |
| `enzymes[].variant_sites` | object | 按变异名映射的候选切点：`{"变异名": [切点...]}`，至少定义 1 个变异名 |

请求示例见 [`examples/panel-request.json`](examples/panel-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/panel \
  -H 'Content-Type: application/json' \
  -d @examples/panel-request.json
```

#### 算法说明（面板）

1. **逐对逐变异评估**：对每一对不同酶、每一个变异型，复用推荐接口的圆环酶切与聚带规则计算参考型 vs 该变异型的条带；任一样本有效切点不足 2 个时该酶对对此变异型不可用。存在至少一条可区分条带时，称该酶对**覆盖**此变异型。
2. **精确最小面板**：把每个酶对的覆盖变异型集合表示为位掩码，按位或 BFS 逐层扩展（第 k 层为恰好 k 个酶对可达的覆盖掩码），首次到达全掩码的层即为最小面板规模——这是最小集合覆盖的**精确解，不使用贪心近似**。同一掩码只保留最优代表酶对（可区分条带总数高者优先，并列取 `pair_key` 小者），同层同掩码只保留最优取值，保证结果完备且确定。
3. **同规模选定**：最小规模相同的方案依次按**可区分条带总数降序**、**排序后的 pair_key 列表字典序**选定。

#### 响应要点（面板）

- `selected_pairs`：最小面板中的酶对（按 `pair_key` 字典序），每个酶对包含
  - `covered_variants`：该酶对覆盖的变异名；
  - `coverages`：**逐变异的差异条带依据**（可区分条带数 + 与推荐接口同结构的 `differences` 明细）；
  - `distinguishable_band_count`：该酶对在全部变异型上的可区分条带总数。
- `covered_variants` / `uncovered_variants`：面板覆盖的变异名 / 无任何可用酶对的变异名。
- `total_distinguishable_bands`：面板全部酶对的可区分条带总数。
- **任一变异型无可用酶对时返回 HTTP 200**：该变异列入 `uncovered_variants`，且 `selected_pairs` 为空——**不返回部分面板**。

### 混合样本解析（POST `/api/v1/digest/mixture`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `components` | array | 命名候选成分列表，**至少 1 个**候选，候选名不可重复，总数不超过 `MAX_MIXTURE_CANDIDATES`（默认 12） |
| `components[].name` | string | 候选成分名（非空白），如 `strain-A` |
| `components[].predicted_bands` | number[] | 该候选成分单独的预测条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `observed_bands` | number[] | 一次凝胶观测到的混合样本条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `band_match_tolerance` | number | ≥ 0，条带容差（bp）：同时用于共迁移聚带与预测/观测配对 |
| `max_components` | int | ≥ 1 且 ≤ 候选总数；枚举成分数 1..该值的全部非空候选子集 |

请求示例见 [`examples/mixture-request.json`](examples/mixture-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/mixture \
  -H 'Content-Type: application/json' \
  -d @examples/mixture-request.json
```

#### 算法说明（混合）

1. **枚举候选子集**：对成分数 1..`max_components` 的全部非空候选子集逐一评估（候选名按字典序规范化，结果与提交顺序无关）。
2. **共迁移聚带（单连接）**：合并子集内各候选的预测条带（多重集），按长度升序扫描，相邻条带长度差 ≤ `band_match_tolerance` 即归入同一簇（链式可传递的单连接聚类）。**共迁移条带只出现一次**，无论它来自多个候选还是同一候选内的重复预测；每簇的**代表长度取簇内全部条带长度的均值**（保留两位小数）。
3. **一对一最小代价匹配**：把合并预测条带与观测条带复用鉴别接口的同一套确定性匹配（先最大化匹配对数，再最小化总误差），得到逐条配对、缺失条带与额外条带。
4. **唯一选定**：全部评估子集依次按 ① **未匹配条带总数**（缺失 + 额外）升序、② **总误差**升序、③ **成分数**升序、④ **排序后的候选名列表字典序**升序比较，次序④保证即使指标完全相同也有唯一结果。

#### 响应示例（混合）

```json
{
  "band_match_tolerance": 10.0,
  "observed_band_count": 4,
  "max_components": 2,
  "selected_components": ["strain-A", "strain-B"],
  "merged_predicted_bands": [230.0, 319.0, 455.0, 700.0],
  "match": {
    "predicted_band_count": 4,
    "matched_count": 4,
    "missing_count": 0,
    "extra_count": 0,
    "total_error": 3.0,
    "matches": [
      {"predicted_length": 230.0, "observed_length": 231.0, "error": 1.0, "description": "..."},
      {"predicted_length": 319.0, "observed_length": 319.0, "error": 0.0, "description": "..."},
      {"predicted_length": 455.0, "observed_length": 455.0, "error": 0.0, "description": "..."},
      {"predicted_length": 700.0, "observed_length": 702.0, "error": 2.0, "description": "..."}
    ],
    "missing_bands": [],
    "extra_bands": []
  },
  "evaluated_subset_count": 3,
  "selection_reason": "在 3 个成分数不超过 2 的候选子集中，选定 ['strain-A', 'strain-B']（2 个成分）……"
}
```

- `selected_components`：选中的候选成分名（字典序）；`merged_predicted_bands`：选中子集聚带后的代表长度（升序）。
- `match`：与鉴别接口**完全同结构**——逐条匹配依据（按预测长度升序）、`missing_bands` / `extra_bands` 及各项计数与总误差。
- `evaluated_subset_count`：实际评估的子集数，等于 `Σ C(候选数, k)`（k = 1..`max_components`）。
- `selection_reason`：中文选定依据，说明四级比较次序如何决定唯一结果。

### 批量盲样分配（POST `/api/v1/digest/batch-identify`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `profiles` | array | 命名预测谱型列表，**至少 1 个**，谱型名不可重复，总数不超过 `MAX_BATCH_PROFILES`（默认 32） |
| `profiles[].name` | string | 谱型名（非空白），如 `strain-A` |
| `profiles[].predicted_bands` | number[] | 该谱型的预测条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `lanes` | array | 命名观测泳道列表，**至少 1 个**，泳道名不可重复，总数不超过 `MAX_BATCH_LANES`（默认 32） |
| `lanes[].name` | string | 泳道名（非空白），如 `lane-1` |
| `lanes[].observed_bands` | number[] | 该泳道的凝胶观测条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `band_match_tolerance` | number | ≥ 0，条带匹配容差（bp）：长度差不超过该值的预测/观测条带才允许配对 |
| `max_unmatched_bands` | int | ≥ 0，单条分配允许的最大未匹配条带数；未匹配数超过该值的“谱型-泳道”边不可用 |

请求示例见 [`examples/batch-identify-request.json`](examples/batch-identify-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/batch-identify \
  -H 'Content-Type: application/json' \
  -d @examples/batch-identify-request.json
```

#### 算法说明（批量盲样分配）

1. **逐边一对一条带匹配**：对每一个“谱型 × 泳道”组合，复用鉴别接口的同一套确定性一对一最小代价匹配（先最大化匹配对数、再最小化总误差），得到该边的**未匹配条带总数**（缺失 + 额外）与**总误差**。未匹配数 **超过** `max_unmatched_bands` 的边标记为不可用，从候选图中剔除（恰好等于阈值的边保留）。
2. **全局一对一最优分配**：在可用边构成的二分图上（每个谱型、每个泳道至多使用一次）用**最小费用最大流**（逐次最短路增广）精确求解，依次：① **已分配泳道数最大化**（最大基数匹配）；② **全部分配边未匹配条带总数最小化**；③ **总误差之和最小化**。
3. **字典序破同分**：三级指标完全相同时，先把泳道按名称排序，再比较该泳道序上“分配到的谱型名序列”的字典序；**未分配项统一排在所有谱型名之后**（名次取谱型总数，比任何真实谱型的 0..P−1 都大）。该字典序被编码为基数 `P + 1` 的位置数并入流的费用（每条分配边对应一个非正增量），因此最优解唯一，且名次仅由名称排序决定，**结果与谱型、泳道的输入提交顺序无关**。

#### 响应示例（批量盲样分配）

```json
{
  "band_match_tolerance": 10.0,
  "max_unmatched_bands": 2,
  "profile_count": 3,
  "lane_count": 3,
  "evaluated_edge_count": 9,
  "eligible_edge_count": 5,
  "assignments": [
    {
      "lane": "lane-1",
      "profile": "alpha",
      "unmatched_count": 0,
      "match": {
        "predicted_band_count": 3,
        "matched_count": 3,
        "missing_count": 0,
        "extra_count": 0,
        "total_error": 7.0,
        "matches": [
          {"predicted_length": 230.0, "observed_length": 231.0, "error": 1.0, "description": "..."},
          {"predicted_length": 320.0, "observed_length": 319.0, "error": 1.0, "description": "..."},
          {"predicted_length": 450.0, "observed_length": 455.0, "error": 5.0, "description": "..."}
        ],
        "missing_bands": [],
        "extra_bands": []
      },
      "description": "泳道 lane-1 分配给谱型 alpha：未匹配条带 0 条（缺失 0 条、额外 0 条），总误差 7 bp……"
    },
    {
      "lane": "lane-3",
      "profile": "beta",
      "unmatched_count": 0,
      "match": { "...": "与鉴别接口完全同结构" },
      "description": "……"
    }
  ],
  "unassigned_lanes": ["lane-2"],
  "unused_profiles": ["gamma"],
  "assigned_lane_count": 2,
  "total_unmatched_count": 0,
  "total_error": 12.0,
  "assignment_reason": "共评估 9 条谱型-泳道候选边……依次以已分配泳道数最大（2 个泳道）、未匹配条带总数最小（0 条）、总误差最小（12 bp）为目标……"
}
```

- `assignments`：逐项匹配依据，**按泳道名字典序排列**；每项含 `lane` / `profile` / `unmatched_count` / `description`，以及与鉴别接口**完全同结构**的 `match`（逐条配对、缺失/额外条带、总误差）。
- `unassigned_lanes`：未分配到任何谱型的泳道名（字典序）；`unused_profiles`：未分配给任何泳道的谱型名（字典序）。
- `evaluated_edge_count` = 谱型数 × 泳道数；`eligible_edge_count` 为未匹配数不超过上限、实际参与全局分配的边数。
- `assigned_lane_count` / `total_unmatched_count` / `total_error`：最优分配的三级目标值；`assignment_reason` 为中文依据，可由响应中的逐项计数复算。

### 自适应决策树（POST `/api/v1/digest/decision-tree`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `plasmid_length` | int | > 0，环状质粒总长度（bp），切点合法范围为 `0 ~ plasmid_length-1` |
| `min_visible_fragment` | int | ≥ 0 且 < 质粒长度；短于此长度的片段不参与条带聚类 |
| `band_resolution_tolerance` | int | ≥ 0，条带分辨容差（bp） |
| `enzymes` | array | 至少 2 种酶，酶名称不可重复；基因型名集合在所有酶间必须完全一致 |
| `enzymes[].name` | string | 酶名称，如 `EcoRI` |
| `enzymes[].genotype_sites` | object | 按基因型名映射的候选切点：`{"基因型名": [切点...]}`，至少定义 1 个基因型 |

基因型数不超过 `MAX_DECISION_GENOTYPES`（默认 12）、酶数不超过独立上限 `MAX_DECISION_ENZYMES`（默认 10）。请求示例见 [`examples/decision-tree-request.json`](examples/decision-tree-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/decision-tree \
  -H 'Content-Type: application/json' \
  -d @examples/decision-tree-request.json
```

#### 算法说明（决策树）

1. **规范化条带签名**：对每一对不同酶（酶名按字典序规范化）× 每一个基因型，复用推荐接口的同一套规则：合并两酶切点（同位点合并）、按圆环相邻切点计算片段（末段 `p0 + L − p_last` 跨原点）、保留长度 ≥ 可见阈值的可见片段，再按分辨容差做**单连接聚带**（相邻片段长度差 ≤ 容差同簇、链式可传递），每簇代表长度取簇内均值（保留两位小数），得到该基因型在该酶对上的规范化签名（代表长度元组，升序）。某基因型在某酶对上有效切点不足 2 个时签名为空元组（与混合接口的共迁移聚带同一纯函数）。
2. **不可区分组合并**：一个基因型在**全部酶对**上的签名构成其签名向量；签名向量完全相同的基因型在任何酶对上都共分支，永远无法区分，先合并为等价类（候选等价类集合）。全部基因型不可区分时树为单个叶节点（`worst_case_depth = 0`）。
3. **记忆化极小化极大搜索**：状态为 `(候选等价类集合, 剩余酶对集合)`。内部节点选择一个**剩余且能拆分当前集合**的酶对做实验，按观测签名把候选类分成若干分支（不能拆分当前集合的酶对直接忽略），递归求解；节点代价为 `1 + max(各分支最坏深度)`。状态做记忆化，因此得到的是**精确的最小最坏实验次数，不是贪心近似**。
4. **同深度破同分**：最坏深度相同时依次取 ① 子树**内部节点总数更少**（一次多分支拆分优先于多层二分支）、② 酶对 **`pair_key` 字典序更小**（酶对按键序枚举，首个最优保留）。基因型划分完全相同的酶对在搜索中完全等价（同为 `available` 时可互换），只保留 `pair_key` 最小的代表参与搜索，不改变最优值。
5. **顺序无关**：全部规范化只依赖切点数据（基因型、酶、切点列表均先排序），分支按签名向量字典序输出，响应与输入提交顺序无关。

#### 响应要点（决策树）

- `genotypes` / `genotype_count`：全部基因型名（字典序）；`total_enzyme_pairs` = C(酶数, 2)。
- `discriminating_pair_count` / `non_discriminating_pairs`：至少能拆分某两个基因型的酶对数，及对全部基因型签名相同、无法拆分任何集合的酶对（按 `pair_key` 字典序）。
- `tree`：决策树根节点（递归结构），每个节点包含
  - `node_id`（按先根遍历从 1 开始）、`depth`（根为 0，每做一次实验 +1）、`node_type`（`internal` / `leaf`）、`candidates`（到达该节点时尚待区分的基因型，字典序）；
  - 内部节点：`pair_key` / `enzyme_pair`（本节点选用的酶对）与 `branches`（按观测签名分组，按签名向量字典序排列；每个分支含 `signature`、落入该签名的 `genotypes` 与递归 `child`）；
  - 叶节点：`indistinguishable` 为 `true` 时该叶是一个永远无法区分的不可区分组（多于 1 个基因型）。
- `worst_case_depth`：最坏实验次数（根到叶最大深度）；`internal_node_count`：实验节点总数；`leaf_count`：叶节点总数。
- `indistinguishable_groups` / `indistinguishable_group_count`：永远无法区分的基因型组（组内与各组均按字典序）；全部基因型可区分时为空。
- `tree_reason`：中文依据，说明签名规范化、不可区分合并、极小化极大目标与两级破同分规则。

### 部分酶切失败反推（POST `/api/v1/digest/partial-failure`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `plasmid_length` | int | > 0，环状质粒总长度（bp）；全部理论片段长度之和必须恰好等于该值 |
| `fragments` | array | 按**圆环顺序**排列的理论片段，至少 1 段，总数不超过 `MAX_PARTIAL_FAILURE_BOUNDARIES`（默认 16） |
| `fragments[].length` | int | 该理论片段长度（bp），> 0 |
| `fragments[].following_boundary` | string | 该片段**后边界**的名称（其后的酶切边界）；各片段后边界名**不可重复**、不可为空白 |
| `observed_bands` | number[] | 一次异常凝胶观测到的条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `band_match_tolerance` | number | ≥ 0，条带匹配容差（bp） |
| `max_failed_boundaries` | int | ≥ 0 且 ≤ 片段总数；枚举失败边界数 0..该值的全部边界子集 |

片段按请求顺序视为圆环顺序：第 i 段 `f_i` 结束于边界 `B_i`，末段 `f_{n-1}` 的后边界 `B_{n-1}` 跨过坐标原点与首段相接。请求示例见 [`examples/partial-failure-request.json`](examples/partial-failure-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/partial-failure \
  -H 'Content-Type: application/json' \
  -d @examples/partial-failure-request.json
```

#### 算法说明（部分失败反推）

1. **枚举失败边界子集**：对失败边界数 k = 0..`max_failed_boundaries` 的**全部** `C(n, k)` 个子集逐一评估（边界按请求中的圆环顺序编号，结果与提交顺序无关）。k = 0 即完全酶切。
2. **失败边界合并**：边界 `B_i` 失败（未切开）意味着片段 `f_i` 末尾不断开、与下一段合并；**连续多个失败边界**把其覆盖的片段并入同一条带；`B_{n-1}` 失败时末段与首段**跨坐标原点**合并。全部边界都失败时整环合并为 1 条长度等于质粒长度的条带。合并后条带长度之和恒等于质粒长度。
3. **一对一条带匹配**：把合并后预测条带（按长度升序）与观测条带复用鉴别接口的同一套确定性一对一最小代价匹配（先最大化匹配对数，再最小化总误差），得到逐条配对、缺失条带与额外条带。
4. **唯一诊断**：全部评估子集依次按 ① **未匹配条带总数**（缺失 + 额外）升序、② **总误差**升序、③ **失败边界数**升序、④ **排序后的失败边界名序列字典序**升序比较，次序④保证即使指标完全相同也有唯一诊断。

#### 响应要点（部分失败反推）

- `failed_boundaries` / `failed_boundary_count`：选定诊断的失败边界名（按名称字典序；完全酶切时为空）及数量。
- `merged_predicted_bands`：合并后预测条带长度（bp），按长度升序。
- `merged_band_details`：合并条带的**组成依据**（按圆环顺序，首条含第 1 个理论片段），每条含
  - `length`：合并条带长度（组成片段长度之和）；
  - `fragment_indexes` / `fragment_lengths`：组成理论片段的序号（按请求顺序从 1 开始）与各自长度；
  - `failed_boundaries_inside`：该条带内部被跨越的失败边界名（条带末尾保留切开的边界不在内），按圆环顺序；
  - `crosses_origin`：该合并条带是否跨过坐标原点（含末段理论片段）。
- `match`：与鉴别接口**完全同结构**——逐条匹配依据（按预测长度升序）、`missing_bands` / `extra_bands` 及各项计数与总误差。
- `evaluated_subset_count`：实际评估的子集数，等于 `Σ C(n, k)`（k = 0..min(上限, n)）。
- `selection_reason`：中文选定依据，说明四级比较次序如何决定唯一结果。

### 凝胶重复实验一致性评估（POST `/api/v1/digest/replicate-consistency`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `lanes` | array | 同一样本的命名观测泳道列表，**至少 1 个**，泳道名不可重复，总数不超过 `MAX_REPLICATE_LANES`（默认 32） |
| `lanes[].name` | string | 泳道名（非空白），如 `lane-1` |
| `lanes[].observed_bands` | number[] | 该泳道的凝胶观测条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `band_match_tolerance` | number | ≥ 0，条带匹配容差（bp）：全局聚类时相邻条带长度差不超过该值即归入同一簇 |
| `min_occurrences` | int | ≥ 1 且 ≤ 泳道数；贡献泳道数达到该值的簇才进入共识谱，否则为低频簇 |

请求示例见 [`examples/replicate-consistency-request.json`](examples/replicate-consistency-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/replicate-consistency \
  -H 'Content-Type: application/json' \
  -d @examples/replicate-consistency-request.json
```

#### 算法说明（重复实验一致性）

1. **全局条带聚类**：把所有泳道的观测条带汇入全局条带池，每条条带记录 `(长度, 泳道序号, 道内序号)`（泳道序号为泳道名字典序名次，道内序号为道内长度升序名次），按该三元组升序排列——**并列按长度与泳道序号破同分**，因此结果与泳道、条带的输入提交顺序无关。扫描有序条带池，相邻条带长度差 ≤ `band_match_tolerance` 即归入同一簇（等价于以容差为阈值的**单连接聚类**，链式可传递）。
2. **共识谱**：每个簇的**出现次数**为贡献该簇的不同泳道数（同一泳道的多条条带只计一次）。出现次数 ≥ `min_occurrences` 的簇进入**共识谱**，共识条带长度取簇内全部条带长度的**中位数**（偶数条取中间两条的均值，保留两位小数）；其余簇为**低频簇**。
3. **逐泳道标记**：对每条泳道——**缺失条带**＝该泳道未贡献条带的共识簇（以共识中位长度记录）；**额外条带**＝该泳道在同一共识簇内除代表条带（与中位长度最接近、并列取较短者）以外的条带；**低频条带**＝该泳道落入低频簇的条带。
4. **稳定性评分**：每条泳道的共识偏差为各命中共识簇 `|代表条带 − 中位长度|` 之和，平均共识偏差＝共识偏差 ÷ 命中共识条带数（未命中时为 0）；异常数＝缺失＋额外＋低频条带数。`稳定性评分 = 100 × 命中共识条带数 ÷（共识条带数 + 额外条带数 + 低频条带数）÷（1 + 平均共识偏差）`，满分 100（无异常且零偏差），异常越多、偏差越大越趋近 0；整体稳定性为各泳道评分的均值。评分只使用响应中输出的计数值与保留两位小数的偏差，可由响应直接复算。

#### 响应示例（重复实验一致性）

```json
{
  "band_match_tolerance": 10.0,
  "min_occurrences": 2,
  "lane_count": 3,
  "total_band_count": 8,
  "cluster_count": 3,
  "consensus_band_count": 3,
  "low_frequency_cluster_count": 0,
  "consensus_profile": [100.0, 200.0, 305.0],
  "clusters": [
    {
      "cluster_id": 1,
      "median_length": 100.0,
      "min_length": 98.0,
      "max_length": 102.0,
      "member_count": 3,
      "occurrence_count": 3,
      "in_consensus": true,
      "lanes": ["lane-1", "lane-2", "lane-3"],
      "members": [
        {"lane": "lane-3", "length": 98.0},
        {"lane": "lane-1", "length": 100.0},
        {"lane": "lane-2", "length": 102.0}
      ]
    }
  ],
  "lanes": [
    {
      "lane": "lane-3",
      "band_count": 2,
      "matched_consensus_count": 2,
      "missing_count": 1,
      "missing_bands": [305.0],
      "extra_count": 0,
      "extra_bands": [],
      "low_frequency_count": 0,
      "low_frequency_bands": [],
      "anomaly_count": 1,
      "consensus_deviation": 12.0,
      "mean_deviation": 6.0,
      "stability_score": 9.52,
      "description": "泳道 lane-3：观测条带 2 条，命中共识条带 2/3，缺失 1 条……"
    }
  ],
  "mean_stability_score": 25.67,
  "score_reason": "全部 8 条观测条带按容差 10 bp 在 3 条泳道间做单连接全局聚类……"
}
```

- `consensus_profile`：共识谱（各共识簇中位长度，升序）；`clusters`：**全部**全局簇明细（含低频簇，`in_consensus` 标记是否入谱），簇编号按中位长度升序。
- `lanes`：逐泳道一致性评估（按泳道名字典序），含缺失 / 额外 / 低频条带明细、异常数、共识偏差与稳定性评分。
- `mean_stability_score`：整体稳定性评分（各泳道评分均值）；`score_reason` 为中文评分依据，可由响应中的计数与偏差复算。

### 拒绝条件（HTTP 422）

以下输入问题整体拒绝，返回统一结构 `{"code":"validation_failed","errors":[...]}`，一次性报告全部错误：

推荐接口（`/api/v1/digest/recommend`）：

- 切点越界（`< 0` 或 `≥ plasmid_length`）；
- 同一样本、同一酶的候选切点重复；
- 酶名称重复；酶少于 2 种；
- 质粒长度 ≤ 0、最小可见片段 ≥ 质粒长度、容差为负等字段约束错误。

鉴别接口（`/api/v1/digest/identify`）：

- 观测条带列表为空（`observed_bands` 至少 1 条）；
- 任一条带长度（参考 / 变异 / 观测）为负数，或匹配容差为负、缺失、非数值；
- 任一条带列表超过 `MAX_BANDS` 上限。

面板接口（`/api/v1/digest/panel`）：

- 参考型或任一变异型的切点越界（`< 0` 或 `≥ plasmid_length`）、同一酶同一样本位点重复；
- 酶名称重复、酶少于 2 种、酶数量超过 `MAX_ENZYMES` 上限；
- 各酶的变异名集合不一致、变异名为空白、未定义任何变异型、变异型数量超过 `MAX_VARIANTS` 上限；
- 质粒长度 ≤ 0、最小可见片段 ≥ 质粒长度、容差为负等字段约束错误。

混合接口（`/api/v1/digest/mixture`）：

- `max_components` 非法（缺失、非整数、< 1）或**超过候选成分总数**；
- 候选成分名**重复**或为空白；候选列表为空；
- 任一候选的 `predicted_bands` 为空（**空候选**），或观测列表为空；
- 候选成分数超过 `MAX_MIXTURE_CANDIDATES`（默认 12）上限；
- 任一条带长度为负数、容差为负 / 缺失 / 非数值等字段约束错误；
- 任一条带列表超过 `MAX_BANDS` 上限。

批量盲样分配接口（`/api/v1/digest/batch-identify`）：

- 谱型或泳道列表为空、谱型名 / 泳道名**重复**或为空白；
- 任一谱型的 `predicted_bands` 或任一泳道的 `observed_bands` 为空；
- 谱型数超过 `MAX_BATCH_PROFILES`（默认 32）或泳道数超过 `MAX_BATCH_LANES`（默认 32）上限；
- `max_unmatched_bands` 缺失、非整数或为负；`band_match_tolerance` 为负 / 缺失 / 非数值；
- 任一条带列表超过 `MAX_BANDS` 上限。

自适应决策树接口（`/api/v1/digest/decision-tree`）：

- 任一基因型的切点越界（`< 0` 或 `≥ plasmid_length`）、同一酶同一基因型位点重复；
- 酶名称重复、酶少于 2 种、酶数量超过独立上限 `MAX_DECISION_ENZYMES`（默认 10）；
- 各酶的基因型名集合不一致、基因型名为空白、未定义任何基因型、基因型数量超过独立上限 `MAX_DECISION_GENOTYPES`（默认 12）；
- 质粒长度 ≤ 0、最小可见片段 ≥ 质粒长度、容差为负等字段约束错误。

部分酶切失败反推接口（`/api/v1/digest/partial-failure`）：

- **理论片段长度之和不等于质粒长度**；
- 任一片段长度 ≤ 0、片段列表为空；
- 后边界名**重复**或为空白；
- 理论片段（边界）数量超过 `MAX_PARTIAL_FAILURE_BOUNDARIES`（默认 16）上限；
- `max_failed_boundaries` 缺失、非整数、为负，或**超过片段（边界）总数**；
- 观测条带列表为空、任一条带长度为负、容差为负 / 缺失 / 非数值；观测列表超过 `MAX_BANDS` 上限。

凝胶重复实验一致性评估接口（`/api/v1/digest/replicate-consistency`）：

- 泳道列表为空、泳道名**重复**或为空白；
- 任一泳道的 `observed_bands` 为空，或任一条带长度为负；
- 泳道数超过 `MAX_REPLICATE_LANES`（默认 32）上限，任一泳道条带数超过 `MAX_BANDS` 上限；
- `band_match_tolerance` 为负 / 缺失 / 非数值；
- `min_occurrences` 缺失、非整数、< 1，或**超过泳道数量**。

注意：推荐接口中切点不足两个属于**酶对级别**的拒绝（记录在 `rejected_pairs`，HTTP 200），而非请求级别错误；面板接口中某变异型无可用酶对时同样返回 HTTP 200（列入 `uncovered_variants`，不返回部分面板）；批量盲样分配接口中所有候选边都超过未匹配阈值时同样返回 HTTP 200（`assignments` 为空、全部泳道列入 `unassigned_lanes`）；决策树接口中某基因型在某酶对上切点不足 2 个或某酶对对全部基因型签名相同，均属于**酶对级别**结果（空签名 / 列入 `non_discriminating_pairs`，HTTP 200），永远无法区分的基因型列入 `indistinguishable_groups`。

---

## 项目结构

```
.
├── app/
│   ├── main.py          # FastAPI 入口：健康检查、推荐/鉴别/面板/混合/批量分配/决策树/部分失败反推/重复实验一致性路由、统一 422 处理
│   ├── models.py        # 请求/响应 Pydantic 模型
│   ├── validation.py    # 越界、重复位点、名称一致性、条带/候选/基因型/泳道规模等语义校验
│   ├── digest.py        # 圆环片段、条带聚类、酶对枚举与排序（纯函数）
│   ├── identify.py      # 条带一对一最小代价匹配与样本判定（纯函数）
│   ├── panel.py         # 逐变异评估与精确最小酶对面板（纯函数）
│   ├── mixture.py       # 混合样本共迁移聚带、子集枚举与最优子集选定（纯函数）
│   ├── batch_identify.py # 批量盲样谱型-泳道边评估与全局一对一最小费用最大流分配（纯函数）
│   ├── decision_tree.py # 规范化条带签名、不可区分合并与记忆化极小化极大决策树（纯函数）
│   ├── partial_failure.py # 失败边界的圆环片段合并（含连续失败/跨原点）、子集枚举与唯一诊断选定（纯函数）
│   ├── replicate.py     # 重复实验全局条带聚类、共识谱、逐泳道缺失/额外/低频标记与稳定性评分（纯函数）
│   └── config.py        # 环境变量配置（均有默认值）
├── tests/               # pytest 单元 + API 端到端测试
├── examples/request.json
├── examples/identify-request.json
├── examples/panel-request.json
├── examples/mixture-request.json
├── examples/batch-identify-request.json
├── examples/decision-tree-request.json
├── examples/partial-failure-request.json
├── examples/replicate-consistency-request.json
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```

## 运行测试

```bash
pip install -r requirements.txt pytest httpx
python3 -m pytest -q
```
