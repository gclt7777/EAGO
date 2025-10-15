# EAGO 算法模块概览

本目录包含基于 PlatEMO 平台实现的 EAGO（Objective Space-based Population Generation to Accelerate Evolutionary Algorithms for Large-scale Many-objective Optimization）算法核心组件。主要文件说明如下：

- `EAGO.m`：算法主入口，负责初始化种群、变量聚类与分析、并在迭代过程中交替执行收敛性优化与分布性优化。文件中包含若干局部函数用于不同阶段的操作（例如 `ConvergenceOptimization_R`、`DistributionOptimization` 等）。
- `EnvironmentalSelection.m`：实现基于非支配排序和角度截断的环境选择，用于在分布性优化后筛选下一代个体。
- `VariableClustering.m`：采用统计量和 K-Means 聚类检测决策变量类型，以区分主导收敛或多样性的变量，为后续优化策略提供输入。

算法通过对变量类型的动态分析来调整不同的优化策略，同时引入多种基于目标空间的生成算子（`GAhalf1/2/3` 等）以提升大规模多目标优化性能。
