# extension-of-HDL-in-genetic-correlation

> ⚠️ 本仓库代码位于 `reference_code/`，是 2020 年前后读博期间的历史实现，仅作参考与溯源，不代表当前可用状态；新的开发请另建目录（如 `src/`），不要直接改这里的文件。

## 这个项目研究什么

这个项目研究如何扩展 HDL 的似然框架，以更精确地估计复杂性状之间的遗传相关与效应量相关。跨 LD 区块回归用于把 SNP 效应汇总到区块层面，检验估计的稳定性。

## 方法/数据

使用 R/Python 编写模拟、跨区块回归及 WLS/GLS 估计脚本，面向 GWAS 汇总统计，并与 HDL、LDSC、GREML 等框架作方法学比较。

## 当前状态

已有跨区块回归和效应量相关的模拟与推导代码；真实数据验证和统一的方法说明尚未完成。

## 代码地图

- `reference_code/cross-block-LD-regression/`：跨 LD 区块回归与相关实现。
- `reference_code/effect-size-correlation/`：效应量相关的模拟、数据整理及不同版本的分析脚本。
