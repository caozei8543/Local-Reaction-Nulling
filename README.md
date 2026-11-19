# Robust MIMO Decoupling via Local Reaction Nulling
## 论文数值复现与真实原始数据分析代码（Python 3.10+）

适用论文：`Robust MIMO Decoupling via Local Reaction Nulling(1).pdf`（5页）及 `RN_Supplementary_Material(1).pdf`（7页）。

本包是根据这两份PDF重新编写的可复现代码，而非作者原始历史代码。它包含可运行的数值复算、数据分析、可视化和约束优化参考实现。**未提供原始HFSS工程、场文件、VNA测量文件、原始优化器名称/起点/历史或具体网格，因此不能宣称独立重跑了文章中的全波仿真或实物实验。** 不会反向制造逐次测量值、误差棒、场分布或仿真曲线。所有 `outputs/reported/` 文件均是PDF已报告数值的转录或代数复算。`tests/` 内的人工测试量仅用于验证算法，不是论文数据。

## 1. 快速运行（无需HFSS）

解压后，在本目录打开终端：

```bash
python -m venv .venv
```

Windows PowerShell：`.venv\Scripts\Activate.ps1`；macOS/Linux：`source .venv/bin/activate`。

```bash
python -m pip install -r requirements.txt
python -m rn reproduce --out outputs/reported
python -m unittest discover -s tests -v
python -m rn --help
```

也可直接执行 `run_reproduction.py`。不依赖 MATLAB、商用软件或联网服务进行上述复算。安装依赖需要正常的 Python 包获取环境。

输出 `Table_I_S1.csv`、`Table_II_S2.csv`、`Table_III_S3.csv`、`Table_S4.csv` 至 `Table_S8.csv`、`summary.json`、Fig.7/S1复数相消重构图和 Fig.10 已报告归一化指标图。后者不包含未知的实测标准差误差棒。

## 2. 核心数值结果（实际运行得到）

| 指标 | 代码结果 | 文章位置 |
|---|---:|---|
| RN归一化反应残差 | 2.003916% | Fig.7 / (S11) |
| 窗内残差比 C_res | 1.391604% | (S12) |
| 反应相位差 | 179.899938° | Fig.7 / (S13) |
| 左/右FD抑制倍数 | 47.38095 / 47.93269 | Table III / S5 |
| 最大复数相对误差 | 5.039029% | Table S5 |
| 最大差模/共模比 | 1.071288% | Table S5 |
| 有限加载数值抑制倍数 | 25.000 / 22.222 | Table III |
| 实测响应抑制倍数 | 20.83333 / 18.18182 | Fig.10 / S7 |
| RN总效率 | 85.769951% | Sec. III-A |
| RN相对带宽 | 6.206897% | Sec. III-A |

复算使用论文中已舍入数字，不能恢复其未舍入精度。对各分母分别处理：反应及FD使用统一复参考 D_ND,L=0.10+j0.05；有限仿真使用同工况数值ND响应；实验使用同工况ND实测平均值。

## 3. 原始数据闭环（需要你已有的真实工程/导出）

将 `inputs/templates/` 拷贝到另一个数据目录，按 `docs/DATA_FORMATS.md` 填入真实数据。模板只有列名/配置说明，没有虚构行。

```bash
python -m rn fields --csv mydata/fields.csv --metadata mydata/field_metadata.json --out results/fields
python -m rn fd --csv mydata/fd.csv --reactions results/fields/field_reactions.csv --out results/fd
python -m rn finite --csv mydata/finite.csv --reactions results/fields/field_reactions.csv --out results/finite
python -m rn vna --csv mydata/vna.csv --out results/vna
python -m rn sweeps --csv mydata/sweeps.csv --out results/sweeps
```

一次运行完整数据分析：

```bash
python -m rn pipeline --manifest mydata/manifest.json --out results/full
```

包含原始数据复算 Fig.3（合适的均匀张量网格时）、Fig.7、Fig.6、Fig.9（提供RN实测扫频时）及 Fig.10（真实样本标准差）。原始 Fig.1、Fig.2、Fig.4、Fig.5、Fig.8 的示意图、CAD、实物照片、HFSS表面电流/场图不是本文程序可以从摘要数据生成的内容，应使用原始工程/图片。

VNA原始文件为标准50 Ω Touchstone v1 `.s2p` 时：

```bash
python -m rn touchstone --manifest mydata/vna_touchstone_manifest.csv --out results/import
python -m rn vna --csv results/import/vna_cycles.csv --out results/vna
```

不插值5.8 GHz，不把 dB 幅度差当作复数差。

## 4. HFSS与结构优化

论文 Stage I 的公共可行初始几何没有给出；`initial_geometry.json` 必须填真实记录。原始优化器算法未公开，本包使用确定性坐标轮询的参考实现，遵守共同可行域、相同初始点、最大25迭代/30目标评估限制。不是原作者的优化轨迹，也不保证找到论文Table I的参数。

`rn/solver.py` 提供完整的ND/DLO/RN目标评估、约束处理、预算记录、固定ND复参考及与真实求解器的JSON进程协议。HFSS工程与求解器端程序须由持有原始工程的人连接。详见 `docs/HFSS_WORKFLOW.md`。对于已有原始HFSS导出，可用 `tools/export_archive_worker.py` 作为精确查询适配器；它只返回已存在的真实导出，任何未存状态直接报错。外部HFSS调用在本环境未实际验证。

```bash
python -m rn optimize --solver mydata/solver_config.json --initial mydata/initial_geometry.json --out results/optimization
python -m rn collect --solver mydata/solver_config.json --out results/experiment
```

`collect` 默认运行Table I三个固定几何，也可以传 `--geometry results/optimization/optimization.json` 验证新优化几何。优化不应在扰动加载工况下重新调谐。

## 5. 文件组织

- `inputs/paper.json`: 来源明确的论文数据和约束；无虚构原始数据。
- `rn/core.py`: 数学核心；复数双线性反应、归一化、差分、相消、效率、可行性。
- `rn/reproduce.py`: 表格和指标的代数重构。
- `rn/raw.py`, `rn/grid.py`: 原始三维场数据、复数FD、有限介质加载和VNA数据分析。
- `rn/sweeps.py`, `rn/touchstone.py`: 频率扫描及标准S2P导入。
- `rn/solver.py`: 约束优化与实际求解器交互协议。
- `rn/plots.py`: 可追溯图表。
- `tests/test_rn.py`: 自动化测试；临时人工夹具不会进入论文输出。
- `docs/TRACEABILITY.md`: 公式、图表及程序对应表与适用边界。
- `docs/DATA_FORMATS.md`: CSV/JSON详细字段规范。
- `docs/HFSS_WORKFLOW.md`: 原始工程验证流程。
- `outputs/reported/`: 当前实际运行产生的复算表图。
- `outputs/TEST_REPORT.txt`: 本次测试日志。

引用本文和原论文既有引用的理论来源。若作为投稿代码附件，请注明它是后续重构分析代码并仅在实际使用了对应数据和软件后报告实际试验结果；不要将参考优化器充当论文历史优化器。投稿作者负责核对原始数据和图表。
