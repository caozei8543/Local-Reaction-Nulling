# 论文—代码—输出追溯表

| 论文内容 | 代码 | 输出 / 边界 |
|---|---|---|
| Eq.(1), S19–S20 | core.material_eta | complex η = εr(1−j tanδ)−1；paper.json |
| Eq.(2), S4–S5 | core.delta_linear, reaction, derivative_from_reaction | 一阶项与真实双线性体积分；S4精确双状态式的计算系数相同，若需直接验证S4需额外导出加载场 |
| Eq.(3), S1–S5 | core.reaction, derivative_from_reaction; raw.fields | 所有三个空间分量及权重；不使用共轭，RMS系数 /2、peak /4 |
| Eq.(4)–(6), S6–S7 | core.common_differential, differential_ratio | 交换对称的数学后果与真实FD的对称性误差 |
| Eq.(7), S11–S13 | core.cancellation; reproduce | 左窗两个复数积分的和、相位差、相消比、幅值失衡贡献 |
| Eq.(8), S7 | core.delta_linear, remainder | 精确D=0时仅在指定幅值子空间中消除线性项；实际RN有残余 |
| S2-D | core.worst_case_first_order | r(|DL|+|DR|)；复数独立幅值的线性上界 |
| Eq.(9), S14–S16 | solver.optimize_one, optimize_all | ND/DLO/RN目标、共同边界及可行域，参考局部算法非历史优化器 |
| Eq.(10), S17–S18 | core.central_difference; raw.fd | 三个步长、复导数、步长减半检查、反应/FD误差 |
| Table I / S1 | inputs/paper.json; reproduce | Table_I_S1.csv（转录） |
| Table II / S2 | reproduce | Table_II_S2.csv（转录） |
| Table III / S3 | reproduce, raw.fd, raw.finite | Table_III_S3.csv（转录＋复相角转换）；真实数据路径另算 |
| Table S4 | reproduce | Table_S4.csv（统一复参考的反归一化） |
| Table S5 | reproduce | Table_S5.csv（复误差与差模/共模比） |
| Table S6 | reproduce | Table_S6.csv（线性预测幅值，不是加载仿真） |
| Table S7 | reproduce | Table_S7.csv（两类有限响应比，分母不同） |
| Table S8 | reproduce | Table_S8.csv（舍入均值代数重构） |
| Fig.3 | raw.fields, plots.height_integrated_density_plot | 真实均匀张量网格时输出高度积分实部；不会伪造场分布 |
| Fig.6 | sweeps.sweeps | 5.5–6.1 GHz真实扫频，采样fmin与最大偏移；不插值未知轨迹 |
| Fig.7 / Fig.S1 | plots.reaction_plot | 已报告复积分重构，或真实场积分版本 |
| Fig.9 | sweeps.sweeps | 原始RN仿真、实测扫频；没有输入不会生成该图 |
| Fig.10 | plots.comparison_plot; raw.vna; pipeline | 报告数据图无未知误差棒；真实数据使用n=10及±1样本SD |
| S21 | raw.vna; touchstone | 每周期先复差后取模，再按工况平均并计算SD |

## 归一化约定

- f0 = 5.8 GHz；e^(+jωt)，无共轭 E1·E2。
- 50 Ω、同参考面、镜像一致端口方向、零相位单位入射功率（不是接受功率）。
- 样本权重 m³、场 V/m；每个子窗54 mm³，每个完整窗108 mm³。
- 对任意 P1/P2，RMS电场系数 −jωε0/(2√(P1P2))；peak 等效系数 −jωε0/(4√(P1P2))。
- 真实流程的归一化基准由ND左窗积分求得，不强行拟合到报告的0.10+j0.05。
- 不将同一窗口内相消解释成两个镜像窗口互相抵消。共同加载时镜像窗口反应相加。
- 论文假设两个固定均匀加载剖面；不能外推到任意位置、任意频率、几何不对称或任意介电对比度。
