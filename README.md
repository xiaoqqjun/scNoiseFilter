# scNoiseFilter

单细胞 RNA-seq 数据噪音基因过滤 R 包

## 概述

scNoiseFilter 是一个用于单细胞 RNA-seq 数据噪音基因识别和过滤的 R 包。支持小鼠、大鼠和人类三种物种，提供 7 类噪音基因的自动识别和过滤功能。

## 安装 / Github

```r
# Install devtools first if needed
install.packages("devtools")

# Install autoSeurat
devtools::install_github("xiaoqqjun/scNoiseFilter")
```

## 快速开始

```r
library(scNoiseFilter)

# 最简单用法 - 直接传入 Seurat 对象
result <- filter_noise_genes(seurat_obj, species = "mouse")

# 查看结果
cat("总基因数:", result$total_genes, "\n")
cat("移除基因:", result$removed_count, "\n")
cat("保留基因:", result$retained_count, "\n")

# 应用过滤
filtered_obj <- subset(seurat_obj, features = result$filtered_genes)
```

## 7 类噪音基因

| # | 类别 | 描述 | 小鼠 | 大鼠 | 人类 |
|---|------|------|------|------|------|
| 1 | **Riken 基因** | 以 `Rik` 结尾的 cDNA 基因 | ✅ | ✅ | ❌ |
| 2 | **预测基因** | Gm 基因 (小鼠) / RGD 基因 (大鼠) | ✅ | ✅ | ❌ |
| 3 | **数字基因** | 含 5+ 连续数字的基因 | ✅ | ✅ | ✅ |
| 4 | **核糖体基因** | Rpl/Rps (小写/大写) | ✅ | ✅ | ✅ |
| 5 | **线粒体基因** | mt-/MT- 前缀 | ✅ | ✅ | ✅ |
| 6 | **血红蛋白基因** | Hba/Hbb/Hbq (小鼠) / HBA/HBB 等 (人类) | ✅ | ✅ | ✅ |
| 7 | **高丰度基因** | 表达丰度中位数超过阈值 | ✅ | ✅ | ✅ |

## 详细用法

### 基本过滤

```r
library(scNoiseFilter)

# 获取基因列表
genes <- rownames(seurat_obj)

# 基本过滤（不包含高丰度过滤）
result <- filter_noise_genes(genes = genes, species = "mouse")
```

### 使用 Seurat 对象（推荐）

```r
# 直接传入 Seurat 对象，自动提取基因和 UMI 矩阵
result <- filter_noise_genes(
  seurat_obj = seurat_obj,
  species = "mouse"
)
```

### 自定义高丰度阈值

```r
# 默认阈值 1% (中位数表达占比)
result <- filter_noise_genes(
  seurat_obj = seurat_obj,
  species = "mouse",
  high_umi_threshold = 0.5  # 0.5%
)
```

### 大数据量采样加速

```r
# 10万+ 细胞时可使用采样加速高丰度计算
result <- filter_noise_genes(
  seurat_obj = seurat_obj,
  species = "mouse",
  sample_rate = 0.1  # 只用 10% 的细胞计算高丰度基因
)
```

### 选择性过滤

```r
# 只过滤特定类型的噪音基因
result <- filter_noise_genes(
  seurat_obj = seurat_obj,
  species = "mouse",
  filter_riken = TRUE,
  filter_predicted = TRUE,
  filter_digits = TRUE,
  filter_ribosomal = FALSE,    # 不过滤核糖体基因
  filter_mito = FALSE,         # 不过滤线粒体基因
  filter_hemoglobin = FALSE,   # 不过滤血红蛋白基因
  filter_high_umi = FALSE      # 不过滤高丰度基因
)
```

### 完整参数示例

```r
result <- filter_noise_genes(
  seurat_obj = seurat_obj,
  species = "mouse",
  sample_rate = 1.0,           # 采样率：1.0=100%, 0.1=10%
  high_umi_threshold = 1,      # 高丰度阈值 (%)
  filter_riken = TRUE,
  filter_predicted = TRUE,
  filter_digits = TRUE,
  filter_ribosomal = TRUE,
  filter_mito = TRUE,
  filter_hemoglobin = TRUE,
  filter_high_umi = TRUE,
  verbose = TRUE
)
```

## 物种支持

### 小鼠 (Mouse)

```r
result <- filter_noise_genes(seurat_obj, species = "mouse")
# 或使用别名
result <- filter_noise_genes(seurat_obj, species = "mmu")
result <- filter_noise_genes(seurat_obj, species = "Mus musculus")
```

支持的噪音基因类别：
- Riken: `Rik[0-9]*$` (如 1700001C02Rik)
- 预测基因: `^Gm\d+$` (如 Gm12345)
- 数字基因: 5+ 连续数字，排除 Riken 和 Gm
- 核糖体: `^Rpl[0-9]|^Rps[0-9]`
- 线粒体: `^mt-|^MT-`
- 血红蛋白: `^Hba|^Hbb|^Hbq`

### 大鼠 (Rat)

```r
result <- filter_noise_genes(seurat_obj, species = "rat")
# 或使用别名
result <- filter_noise_genes(seurat_obj, species = "rno")
```

支持的噪音基因类别：
- Riken: `Rik[0-9]*$`
- 预测基因: `^RGD\d+$` (如 RGD12345)
- 其他与小鼠相同

### 人类 (Human)

```r
result <- filter_noise_genes(seurat_obj, species = "human")
# 或使用别名
result <- filter_noise_genes(seurat_obj, species = "hsa")
```

支持的噪音基因类别（无 Riken 和预测基因）：
- 数字基因: 5+ 连续数字
- 核糖体: `^RPL[0-9]|^RPS[0-9]`
- 线粒体: `^MT-|^MTND|^MTCO|^MTATP|^MTCYB`
- 血红蛋白: `^HBA|^HBB|^HBQ|^HBD|^HBE|^HBG|^HBM|^HBZ`

## 返回值

`filter_noise_genes()` 返回一个列表：

```r
result <- filter_noise_genes(seurat_obj, species = "mouse")

# 返回值结构
result$total_genes      # 输入基因总数
result$noise_genes     # 各类别噪音基因详情
result$filtered_genes  # 过滤后保留的基因
result$removed_genes   # 所有被移除的基因
result$removed_count   # 移除基因数量
result$retained_count  # 保留基因数量
result$species         # 使用的物种
result$params          # 使用的参数

# 查看各类别详情
result$noise_genes$riken      # Riken 基因
result$noise_genes$predicted  # 预测基因
result$noise_genes$digits     # 数字基因
result$noise_genes$ribosomal  # 核糖体基因
result$noise_genes$mito       # 线粒体基因
result$noise_genes$hemoglobin # 血红蛋白基因
result$noise_genes$high_umi   # 高丰度基因
```

### 各类别详情结构

```r
# 每个类别包含
result$noise_genes$riken$name      # 类别名称
result$noise_genes$riken$pattern   # 使用的正则表达式
result$noise_genes$riken$genes     # 匹配的基因列表
result$noise_genes$riken$count     # 匹配的基因数量

# 高丰度基因额外包含 details
result$noise_genes$high_umi$details  # data.frame: gene, median_percent
```

## 辅助函数

### get_umi_matrix()

从 Seurat 对象提取 UMI 矩阵，自动兼容 v4 和 v5：

```r
# 自动检测 Seurat 版本并提取
umi_matrix <- get_umi_matrix(seurat_obj)
```

### get_noise_patterns()

获取物种特异的噪音基因 pattern：

```r
patterns <- get_noise_patterns("mouse")
patterns$riken$pattern     # "Rik[0-9]*$"
patterns$predicted$pattern # "^Gm\\d+$"
```

### get_noise_report()

生成过滤报告：

```r
report <- get_noise_report(result, format = "list")
# format: "text", "data.frame", "list"
```

## 高丰度基因计算方法

高丰度基因使用**表达丰度中位数法**：

1. 计算每个基因在每个细胞中的 UMI 占比（百分比）
2. 对每个基因，取所有细胞的中位数
3. 中位数超过阈值的基因被标记为高丰度

```r
# 阈值含义示例
# threshold = 1: 基因在超过一半的细胞中，每个细胞至少 1% 的 UMI
# threshold = 0.5: 基因在超过一半的细胞中，每个细胞至少 0.5% 的 UMI

# 推荐阈值
# 2%: 极严格，过滤极高表达基因（如严重背景污染）
# 1%: 严格，过滤高表达看家基因（如 Malat1）
# 0.5%: 默认，过滤较高表达基因（如 Actb）
```

## 采样加速（sample_rate）

大数据量时，可使用采样加速高丰度基因的计算：

### 参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `sample_rate` | `1.0` | 采样比例，1.0 = 100%（不采样） |

### 推荐值

| 细胞数 | 推荐 sample_rate | 说明 |
|--------|-----------------|------|
| < 5万 | `1.0` (100%) | 不需要采样 |
| 5-10万 | `0.5` (50%) | 采样一半 |
| 10-50万 | `0.2` (20%) | 采样 20% |
| > 50万 | `0.1` (10%) | 采样 10% |

### 使用示例

```r
# 小数据量 - 不采样（默认）
result <- filter_noise_genes(seurat_obj, species = "mouse")

# 大数据量 - 采样 10% 细胞加速
result <- filter_noise_genes(
  seurat_obj = seurat_obj,
  species = "mouse",
  sample_rate = 0.1  # 10% 采样
)
```

### 注意事项

- 采样**仅影响高丰度基因计算**，不影响其他 6 类 pattern 过滤
- 采样使用固定随机种子 (`set.seed(42)`)，结果可重复
- 采样后结果与全量计算差异很小（高丰度基因通常稳定）

## 与 Seurat 兼容性

scNoiseFilter 兼容 Seurat v4 和 v5：

```r
# 自动检测版本并提取数据
umi_matrix <- get_umi_matrix(seurat_obj)

# 内部实现
# Seurat v5: Seurat::LayerData(obj, layer = "counts")
# Seurat v4: Seurat::GetAssayData(obj, slot = "counts")
```

## 完整工作流示例

```r
library(Seurat)
library(scNoiseFilter)

# 加载 Seurat 对象
seurat_obj <- readRDS("my_seurat_object.rds")

# 1. 运行噪音基因过滤
result <- filter_noise_genes(
  seurat_obj = seurat_obj,
  species = "mouse",
  high_umi_threshold = 0.5,
  verbose = TRUE
)

# 2. 查看结果摘要
cat("原始基因数:", result$total_genes, "\n")
cat("移除基因数:", result$removed_count, "\n")
cat("保留基因数:", result$retained_count, "\n")
cat("过滤率:", round(result$removed_count / result$total_genes * 100, 1), "%\n")

# 3. 查看各类别
for (category in names(result$noise_genes)) {
  info <- result$noise_genes[[category]]
  cat(sprintf("  %s: %d genes\n", info$name, info$count))
}

# 4. 应用过滤
filtered_obj <- subset(seurat_obj, features = result$filtered_genes)

# 5. 保存结果
saveRDS(filtered_obj, "filtered_seurat_object.rds")
```

## 注意事项

1. **大小写敏感**: 基因名匹配区分大小写（如 Grik1 不会被当作 Riken 基因）
2. **高丰度过滤**: 需要传入 Seurat 对象或 UMI 矩阵才能启用
3. **物种选择**: 人类没有 Riken 和预测基因类别

## 版本历史

### v1.0.0
- 初始版本
- 支持小鼠、大鼠、人类三种物种
- 7 类噪音基因过滤
- Seurat v4/v5 兼容
- 高丰度基因中位数法
- 采样加速支持 (sample_rate)

## 作者

**Zhijun Feng (冯志军)**

- **邮箱**: fenzhj18@sina.com / xiaoqqjun@sina.com
- **GitHub**: https://github.com/xiaoqqjun/scNoiseFilter
- **ORCID**: 0000-0003-1813-1669
- **微信**: 博士后的小酒馆

## 许可证

MIT License
