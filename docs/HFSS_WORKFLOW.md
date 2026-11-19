# HFSS 原始工程与自动化求解接口

## 1. 为什么没有 .aedt 工程

论文提供了基板 RO4350B、εr=3.66、tanδ=0.0037、厚0.762 mm、40×20 mm²、两个9×6×2 mm³窗、50 Ω等信息及Table I几何参数。但不提供完整原始CAD实体、馈线/连接器模型、窗绝对坐标、边界网格、求解设置、优化器名称与起点。没有这些信息不能制造与真实文章完全相同的HFSS工程。请使用作者存档工程。

## 2. 建议实际HFSS工作流

1. 打开已建立的原始镜像对称双端口模型；核对 Fig.2/Table I 的变量含义（mm）。使用同一基板、材料、地面、两端口参考平面和归一化。设置实际两个加载窗相对于空气的介电参数。不要误将窗内两个反应子区赋为不同介质。
2. 分别输入ND、DLO、RN Table I几何，名义空气状态求解5.8 GHz。检查 |S21|≤−20 dB、两端口 |Sii|≤−15 dB、辐射效率≥0.8。导出S参数复实虚及辐射效率；保存HFSS版本与设置。
3. 对每个名义解分别激励Port 1、Port 2，各1 W入射功率、0°，另一端口匹配。导出同坐标三个复电场分量，确认HFSS的peak/RMS约定。使用 `grid` / `merge-fields` / `fields` 计算反应和D，并做网格细化。
4. 对每设计在左/右窗分别设置 ±0.01、±0.005、±0.0025 实介电扰动，不改变几何及频率。保存全部复S21。运行 `fd` 验证反应预测和对称性。
5. 按εr=ε′(1−j tanδ)的损耗约定加载51/51、51/71 HF，保持窗位置和原始几何，导出复S与效率。运行 `finite` 计算高阶余量和有限加载响应。对左窗另执行−0.1至0.1扫频，运行 `sweeps`。
6. 用原始VNA记录运行 `vna`。真实空载变化记录独立保存，不将理论或仿真当作测量。

## 3. 自动化接口（任意真实HFSS控制脚本）

`rn/solver.py` 不依赖某一版本的 PyAEDT 或 COM 接口，它通过已经实现的本地进程协议调用用户实际的HFSS worker。这样避免未验证API造成错误，也不需要在没有完整工程时猜测端口对象名。将 `inputs/templates/solver_config.json` 中 command 替换为真实脚本的 argv 列表（不使用 shell）。worker 接收两个参数：请求JSON的绝对路径、响应JSON的绝对路径。`solver.py` 会等待完成后验证返回请求与原请求完全相等。

请求字段（所有数值都由程序填好）：

```json
{"geometry_mm":{"l1":1.5,"l2":1.5,"d1":1.4,"d2":1.4,"ws":1.65,"xf":6.0,"lp":11.0},"eta_l_re_im":[0.0,0.0],"eta_r_re_im":[0.0,0.0],"frequency_hz":5800000000.0,"reference_ohm":50,"need_reaction":true}
```

Worker责任：以实际HFSS工程设置请求几何和负载，求解5.8 GHz。atan等不得将介质损耗号搞反。请求η为复数；材料εr=1+η，若以实部和tanδ输入，可计算 `eps_real=1+eta.real`、`tan_delta=-eta.imag/eps_real`。其他设置不变。`need_reaction=true` 时必须导出两个名义端口激励场并做空间积分，返回 **RMS且除以√(P1P2)的SI复反应**，不可返回归一化至ND后的相对值。若原始场为peak，则先将两个场各除√2。端口参考一致。

输出JSON字段（括号描述，不是数据示例）：

- `status`: `ok`
- `request`: 原请求完整JSON（不能改变数值）
- `provenance`: `full_wave_export`（仅当确为真实全波结果）
- `s11_re_im`, `s21_re_im`, `s22_re_im`: [实部,虚部]；均为50 Ω线性复散射参数
- `radiation_efficiency`: 0..1
- `reaction_l_SI_re_im`, `reaction_r_SI_re_im`: [实部,虚部]；请求需要时必填

worker发生错误必须返回非零退出状态或不含ok的响应。不要填充默认仿真值。`ExternalSolver` 的 `worker_requests` 是实际请求次数，而非内部HFSS求解次数或耗时。可选的 `tools/export_archive_worker.py` 使用已有真实导出索引（每条完整response），精确匹配请求并返回响应，支持已有状态的复算；不能求解新几何，因此不能替代完整HFSS优化。配置 command 可设为 `["python", "tools/export_archive_worker.py", "mydata/real_export_index.json"]`。真实索引为JSON数组，每项上述完整响应。

## 4. 优化说明

`optimize` 使用共同、已经可行的起点（需要原记录），几何边界来自 Table I。ND目标 |S21(0)|²；DLO为左右窗各±0.05/±0.10的复差模平方和；RN目标 |RL/RL_ND|²。先得到ND最优可行结果并单独查询固定参考，再对同一起点执行DLO、RN。每设计最多25次轮询迭代、30次候选目标评估；重复几何缓存；不合格候选计入预算却不接受。相同目标评估次数并不等价于相同全波求解次数。程序日志清楚区分这些计数。不同局部算法不一定重现Table I的尺寸，故原论文历史结果需要原始优化器记录独立验证。

真实HFSS自动化worker因依赖原始工程以及特定软件安装，未在本环境执行；未声称已完成全波求解测试。
