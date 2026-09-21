"""运行期配置：全部通过环境变量注入，均有默认值，核心功能不依赖任何外部服务。"""
import os


class Settings:
    APP_NAME: str = os.getenv("APP_NAME", "circular-dna-digest")
    # 服务监听地址与端口（README / Compose 中使用 8000）
    HOST: str = os.getenv("APP_HOST", "0.0.0.0")
    PORT: int = int(os.getenv("APP_PORT", "8000"))

    # 单次请求允许提交的最大酶数量，避免组合枚举过大
    MAX_ENZYMES: int = int(os.getenv("MAX_ENZYMES", "64"))

    # identify 接口每个条带列表允许的最大条带数，控制匹配 DP 的规模
    MAX_BANDS: int = int(os.getenv("MAX_BANDS", "512"))

    # panel 接口允许的最大变异型数量：最小面板为精确集合覆盖，规模随变异数指数增长
    MAX_VARIANTS: int = int(os.getenv("MAX_VARIANTS", "16"))

    # mixture 接口允许的最大候选成分数：需枚举至多该规模的全部子集
    MAX_MIXTURE_CANDIDATES: int = int(os.getenv("MAX_MIXTURE_CANDIDATES", "12"))

    # batch-identify 接口允许的最大谱型数与泳道数：候选边数为两者乘积，
    # 全局一对一分配为最小费用最大流
    MAX_BATCH_PROFILES: int = int(os.getenv("MAX_BATCH_PROFILES", "32"))
    MAX_BATCH_LANES: int = int(os.getenv("MAX_BATCH_LANES", "32"))

    # decision-tree 接口独立的规模上限：极小化极大搜索的记忆化状态为
    # “候选基因型集合 × 剩余酶对集合”，规模随基因型数与酶对数快速增长
    MAX_DECISION_GENOTYPES: int = int(os.getenv("MAX_DECISION_GENOTYPES", "12"))
    MAX_DECISION_ENZYMES: int = int(os.getenv("MAX_DECISION_ENZYMES", "10"))

    # partial-failure 接口允许的最大理论片段（边界）数：需枚举至多
    # max_failed_boundaries 规模的全部边界子集（最坏 2^n），随边界数指数增长
    MAX_PARTIAL_FAILURE_BOUNDARIES: int = int(
        os.getenv("MAX_PARTIAL_FAILURE_BOUNDARIES", "16")
    )

    # replicate-consistency 接口允许的最大观测泳道数：全局聚类规模为
    # 泳道数 × 每泳道条带数（每条带列表仍受 MAX_BANDS 约束）
    MAX_REPLICATE_LANES: int = int(os.getenv("MAX_REPLICATE_LANES", "32"))


settings = Settings()
