# 原始数据输入规范

所有 CSV 均 UTF-8，首行为列名。实部/虚部对应复数线性量，不能填 dB 或角度。输入必须来自原始软件导出而非拟合。所有固定频率记录必须有恰好5.8 GHz样点（容差1 Hz）。字段模板在 `inputs/templates`。`design` 取 ND / DLO / RN。文件中数值用小数点。

## 场积分

`fields.csv`: design, frequency_hz, region (La/Lb/Ra/Rb), x_m,y_m,z_m, weight_m3,g, E1x_re,E1x_im,...,E2z_re,E2z_im。每行是同一空间坐标下两个独立端口激励的复数电场。区域不得交叠，四个子区域各54 mm³。均匀profile取g=1。每个design/region必须完整；权重不重复，代码检查各区域总体积误差≤1%。

`field_metadata.json`: 确认可比性而非猜测，包含时间约定、端口50 Ω、零相位、另一端口匹配、单位、phasor `rms` / `peak`、两个入射功率。模板的peak必须由实际软件确认。若是 accepted-power 数据，不得直接导入。两次激励必须在相同电磁模型、网格、坐标、端口参考面完成。坐标/权重来源和场导出原始文件另保存。

用 `grid` 从真实CAD坐标生成中点等体积积分网格；四个区域定义必须由实际模型测得（不得按示意图猜位置）。用 `merge-fields` 将同一网格两个端口场文件对齐：

```bash
python -m rn grid --config mydata/grid_config.json --out results/grid
python -m rn merge-fields --grid results/grid/quadrature_grid.csv --port1 ND_p1.csv --port2 ND_p2.csv --design ND --out results/ND
```

端口场CSV列名 x_m,y_m,z_m,Ex_re,Ex_im,Ey_re,Ey_im,Ez_re,Ez_im；坐标导出时保持相同文本精度，代码精确匹配坐标。对DLO、RN重复，然后将三个生成的fields.csv仅保留一行表头纵向合并。使用至少两套网格复算RN残差，比较归一化之前的绝对复积分，报告幅值和相位变化；PDF仅报道汇总值，代码不能恢复真实收敛历史。Fig.3式图在均匀张量网格下沿z(高度)积分；若真实网格不规则，应先用实际体积权重/几何构造二维投影，不应伪造像素。

## 复数FD

`fd.csv`: design,frequency_hz,window,delta,sign,s21_re,s21_im。window=L/R；delta=0.01,0.005,0.0025；sign=+/-。3×2×3×2=36个状态。每次仅改变指定窗口的 εr=1±delta（另一窗εr=1），保持频率和几何不变。负扰动是数值测试，不是泡沫实验。读取的复场/复S应使用同一50 Ω参考。步长变化以细步长幅值作分母，使用主值相差；原文未给出更细的计算约定，需记录此定义。

## 有限介质加载

`finite.csv`: design,condition,frequency_hz,s21_re,s21_im,s11_re,s11_im,s22_re,s22_im,radiation_efficiency,eta_l_re,eta_l_im,eta_r_re,eta_r_im。condition=baseline/common/unequal。baseline=空气；common=51/51 HF；unequal=51/71 HF（L为51、R为71的约定；原论文未指定不等载荷具体左右顺序，理论镜像交换等价）。每种9个总状态（3设计×3工况）。εr复扰动以原文5 GHz代表料参用于5.8 GHz。efficiency=0..1小数。软件输出可能用不同的accepted-power效率定义，须核对。比较同时保留|ΔS21|、高阶余量复数与模、匹配与效率。程序不会根据两个归一化幅值逆推高阶余量相位。

## VNA与误差棒

`vna.csv`: design,condition,cycle,frequency_hz,loaded_re,loaded_im,unloaded_re,unloaded_im。condition=common/unequal；cycle=1..10，每工况每设计10行，总计60行。unloaded 必须该周期实际使用的复数参考，允许同一真实公用参考重复填入。真实参考方式须写入实验记录。标准差为样本SD (ddof=1)，非标准误差。`--unloaded-bound` 仅在原始空载变化记录支持所填数值时使用。

`vna_touchstone_manifest.csv`: design,condition,cycle,loaded_s2p,unloaded_s2p；文件路径相对于manifest。须为Touchstone v1单端50 Ω S参数，支持 RI/MA/DB；读取顺序S11,S21,S12,S22；不支持v2或端口混合模；必须存在5.8 GHz样点。测量导出应保留足够小数位以解析RN约10^-4量级的变化。

## 频率扫描

`sweeps.csv`: design,source,eta_l,eta_r,frequency_hz,s11_re,s11_im,s21_re,s21_im。source=simulation/measurement；这里只允许实数η的频率图。Fig.6 要求三设计左窗扰动−0.1至0.1内、eta_r=0且包含eta_l=0的扫频，共同频率网格；使用原始采样最小值不做二次插值。Fig.9 要另有RN measurement、eta_l=eta_r=0的标称扫频。输出不补齐缺失载荷。

## 数据留存

保留实际HFSS版本、EDA工程哈希、馈线/SMA参考校准、网格及端口设置、求解频点、RF源位相、VNA导出文件、实际加载位置及顺序、优化器配置/终止历史。原文未给出的信息不能由代码生成。若输入数据不完整，对应功能明确报错。
