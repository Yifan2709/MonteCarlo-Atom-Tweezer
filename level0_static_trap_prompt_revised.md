# Codex 执行 Prompt：Level 0 数值框架与单光镊径向静态检查（修订版）

请在当前工作区中实际完成“Level 0：数值框架与单光镊径向静态检查”。不要只给方案或代码片段；必须创建结构化代码、安装或确认依赖、运行测试和完整模拟，并生成数据、图表与汇总报告。

如果工作区存在 `sources/`，将其中所有内容视为只读参考材料，不得编辑、移动、重命名或删除，并在结束前检查其工作区状态；如果 `sources/` 不存在，直接继续，不要为满足本条而创建空目录。

## 一、任务定位与边界

本任务受 arXiv:2403.12021v4 的实验参数启发，但不是论文复现。它只建立一个静态、一维、经典、单原子的径向光镊模型，用作后续 SLM-AOD 转移和运输模拟之前的数值自洽检查。

必须在 README、`metrics.json` 和 `summary.md` 中明确以下限制：

- (x) 是高斯光束焦平面内的一个径向坐标，不是光轴方向；
- 不模拟论文报告的轴向 5.64 kHz 运动；
- 不模拟 SLM-AOD 交接、深度随时间变化、长距离移动、AOD 柱面透镜效应、二维运动、噪声、量子内态、加热、损失或多原子；
- SLM 与 AOD 在本 Level 0 中调用同一高斯势函数，二者只通过配置参数和参数来源区分；
- 数值曲率与解析曲率的比较是代码一致性检查，不是独立实验验证。

分别完成以下内容：

1. 建立一维径向高斯光偶极势；
2. 计算对应受力；
3. 生成 SLM 和 AOD 的静态势阱；
4. 数值寻找势阱平衡位置并检查力残差；
5. 计算势阱中心曲率和小振幅频率；
6. 比较解析曲率、数值曲率和相应频率；
7. 使用 velocity-Verlet 模拟静态势阱中的一条经典轨迹；
8. 计算动能、势能、总能量和归一化能量误差；
9. 从轨迹估计有限振幅振动频率；
10. 将离散化误差与有限振幅非简谐频移明确分开；
11. 生成 CSV、JSON、Markdown 和 PNG 结果；
12. 通过自动化测试与程序化图像检查验收结果。

## 二、物理模型

在光束焦平面的一个径向方向上使用

\[
U(x)=-U_0\exp\!\left[-\frac{2(x-x_c)^2}{w^2}\right],
\]

其中 (U_0>0) 是势阱深度，(x_c) 是势阱中心，(w) 是标准的 (1/e^2) 强度半径。

力为

\[
F(x)=-\frac{\mathrm dU}{\mathrm dx}
=-\frac{4U_0(x-x_c)}{w^2}
\exp\!\left[-\frac{2(x-x_c)^2}{w^2}\right].
\]

中心解析曲率、小振幅角频率和普通频率分别为

\[
\kappa_0=\left.\frac{\mathrm d^2U}{\mathrm dx^2}\right|_{x=x_c}
=\frac{4U_0}{w^2},
\qquad
\omega_0=\sqrt{\frac{\kappa_0}{m}},
\qquad
f_0=\frac{\omega_0}{2\pi}.
\]

经典运动和能量为

\[
m\ddot{x}=F(x),
\qquad
E(t)=\frac12mv(t)^2+U[x(t)].
\]

势能零点取在无穷远处。只有 (E_0<0) 的轨迹是束缚轨迹；运行前必须检查这一条件。如果初始条件不束缚，给出清晰中文错误并停止该势阱的轨迹频率分析。

### 参数与论文的关系

- SLM 默认深度 180 μK 对应论文报告的阵列满深度 (k_B\times0.18(2)\) mK；
- AOD 默认深度 280 μK 对应论文运输/转移使用的 AOD 深度；
- `waist_um: 1.17` 是论文对 1061 nm SLM 阵列报告的束腰；
- 论文没有单独报告 AOD 束腰。AOD 默认也取 1.17 μm 只是因为其后孔径光束尺寸被匹配到 SLM 光路，必须在配置、报告和图表元数据中标记为 `assumed`；
- 180 μK SLM 与 280 μK AOD 是两个独立静态势阱的 Level 0 比较，不代表论文某一时刻的实际重叠转移势。论文 Fig. 6 转移时的 SLM 深度约为 140 μK；
- 论文对 1061 nm SLM 测得的径向和轴向频率分别是 29.30(4) kHz 和 5.64(3) kHz。Level 0 只把径向值作为量级参照，不把轴向值作为当前模型的验收目标。

## 三、单位与默认配置

程序内部统一使用 SI 单位。集中管理以下换算：

\[
U[\mathrm J]=k_B U[\mu\mathrm K]\times10^{-6},
\quad
w[\mathrm m]=w[\mu\mathrm m]\times10^{-6},
\quad
m[\mathrm{kg}]=m[\mathrm u]m_u.
\]

默认 YAML 使用以下结构：

```yaml
model:
  coordinate: radial_at_focus

atom:
  isotope: Cs-133
  mass_u: 132.90545196

slm:
  name: SLM
  depth_uK: 180.0
  waist_um: 1.17
  center_um: 0.0
  depth_source: paper_full_depth
  waist_source: paper_1061nm_slm_measurement

aod:
  name: AOD
  depth_uK: 280.0
  waist_um: 1.17
  center_um: 0.0
  depth_source: paper_transport_depth
  waist_source: assumed_matched_to_slm_back_aperture

grid:
  min_offset_waist: -4.0
  max_offset_waist: 4.0
  num_points: 4001
  curvature_step_waist: 0.002

trajectory:
  initial_displacement_waist: 0.10
  initial_velocity_m_per_s: 0.0
  duration_us: 3000.0
  dt_us: 0.05
  frequency_estimator: interpolated_zero_crossings
  minimum_complete_cycles: 20

validation:
  max_omega_dt: 0.05
  curvature_relative_tolerance: 1.0e-4
  harmonic_frequency_relative_tolerance: 1.0e-4
  trajectory_frequency_relative_tolerance: 0.03
  max_energy_error_over_depth: 1.0e-4

output:
  directory: outputs/level0_static_trap/demo
  save_png: true
  save_csv: true
  save_json: true
```

说明：

- `initial_displacement_waist` 是相对势阱中心的位移，必须按 (x(0)=x_c+q_0w) 解释；不得把它当成实验室坐标；
- 网格上下限也是相对中心的无量纲偏移；
- `duration_us: 3000` 给出约 87 个 SLM 周期，原始 FFT 栅格约 0.333 kHz。默认频率估计仍使用插值过零点，FFT 只作诊断；
- 若用户将时长改短，必须根据解析频率检查完整周期数；少于 `minimum_complete_cycles` 时不得声称获得了高精度轨迹频率；
- `trajectory_frequency_relative_tolerance` 不是积分器误差容差。有限振幅高斯势具有非简谐软化，轨迹频率相对小振幅频率出现约 1% 量级、在更大允许振幅下达到 1-3% 的差异是物理效应；
- 若要把轨迹频率用于严格的谐振子回归测试，测试配置应使用 `initial_displacement_waist: 0.01`，不得用生产默认值 0.10 强行要求 (10^{-4}) 级一致。

配置加载必须检查所有数值是否有限，并验证：质量、深度、束腰、时长、步长均大于零；网格点数至少为 5；偏移范围有序；曲率为正；时间数组严格递增且等间隔；初始状态束缚；预估完整周期数足够；

\[
\omega_{\max}\Delta t < \texttt{max\_omega\_dt}.
\]

若不满足，给出包含实际值和阈值的中文错误信息。

## 四、项目结构与安装方式

使用可安装的 `src` 布局，必须包含 `pyproject.toml`，不能只创建 `requirements.txt` 后假定包可导入：

```text
simulation/
└── level0_static_trap/
    ├── README.md
    ├── pyproject.toml
    ├── requirements.txt
    ├── configs/
    │   └── default.yaml
    ├── src/
    │   └── level0_static_trap/
    │       ├── __init__.py
    │       ├── constants.py
    │       ├── config.py
    │       ├── potentials.py
    │       ├── analysis.py
    │       ├── integrators.py
    │       ├── simulation.py
    │       ├── visualization.py
    │       ├── io_utils.py
    │       └── cli.py
    ├── tests/
    │   ├── test_config.py
    │   ├── test_potentials.py
    │   ├── test_frequency.py
    │   ├── test_integrator.py
    │   └── test_cli_smoke.py
    └── outputs/
        └── .gitkeep
```

`pyproject.toml` 至少应：

- 使用 setuptools 或其他轻量标准构建后端；
- 声明 Python 最低版本；
- 声明 NumPy、SciPy、Matplotlib 和 PyYAML 运行依赖；
- 将 pytest 声明为开发依赖；
- 正确配置 `src` 包发现；
- 可选地提供 `level0-static-trap` console script。

从 `simulation/level0_static_trap/` 目录执行并实际验证：

```bash
python -m pip install -e .
python -m pytest
python -m level0_static_trap.cli \
  --config configs/default.yaml \
  --trap both \
  --output-dir outputs/level0_static_trap/demo
```

最终回复不得声称 CLI 可运行，除非上述模块命令确实以退出码 0 执行。不得用临时修改 `PYTHONPATH` 掩盖打包缺失。

## 五、模块职责

### `constants.py`

集中定义或从权威常数库集中导入：(k_B)、(m_u) 和必要单位换算。其他模块不得重复定义。

### `config.py`

- 读取和验证 YAML；
- 拒绝未知的关键字段或至少明确警告；
- 转换为结构化配置对象和 SI 单位；
- 保留输入值、单位、来源及假设标签；
- 解析输出路径时给出明确规则，不能依赖任意当前工作目录；
- CLI 的 `--output-dir` 覆盖值必须写入 `config_used.yaml`。

### `potentials.py`

至少实现：

```python
gaussian_potential(x, depth, waist, center)
gaussian_force(x, depth, waist, center)
gaussian_curvature_analytic(depth, waist)
trap_frequency_analytic(depth, waist, mass)
```

支持标量和 NumPy 数组；所有公开函数具有简短中文 docstring；SLM 和 AOD 不得复制实现。

### `analysis.py`

至少实现：

```python
find_equilibrium_position(...)
numerical_curvature(...)
numerical_trap_frequency(...)
estimate_frequency_from_trajectory(...)
estimate_nonlinear_frequency_by_quadrature(...)
calculate_relative_error(...)
summarize_static_trap(...)
```

要求：

1. 平衡位置搜索使用有界、无量纲化的坐标或等价稳定方法；不得直接依赖焦耳量级目标函数的默认绝对容差；
2. 返回平衡位置、势能和力残差；
3. 曲率至少使用 (h) 和 (h/2) 两个步长，报告收敛差异；
4. 曲率频率与解析频率的误差归类为数值离散化误差；
5. 轨迹频率默认由同方向、带插值的过零点或多峰线性回归估计；
6. FFT 可以输出诊断频谱，但必须同时报告 (1/T) 的原始频率分辨率，不得仅按最高 FFT bin 声称百分之一精度；
7. 通过一维能量积分计算同一初始能量下的非线性周期，或等价地明确报告有限振幅频移；
8. 分别输出：解析小振幅频率、数值曲率频率、非线性积分基准频率和轨迹估计频率；
9. `calculate_relative_error` 的参考量为零时必须显式处理。

### `integrators.py`

实现针对静态、位置依赖力的 velocity-Verlet：

```python
velocity_verlet(force_function, mass, x_initial, v_initial, times)
```

返回时间、位置、速度和加速度数组。必须检查时间数组等间隔。

本 Level 0 只保证自治势 (F(x)) 下的辛性质。README 应说明：未来若扩展到 (F(x,t))，接口、能量解释和辛结构都需要重新审视；不得仅凭“传入函数”就宣称当前积分器已经正确支持任意时变势。

### `simulation.py`

对每个选中的势阱使用同一流程：

1. 建立网格；
2. 计算势能和力；
3. 求平衡位置和力残差；
4. 计算解析和数值曲率；
5. 计算解析、数值曲率及非线性基准频率；
6. 检查时间步长、周期数和束缚条件；
7. 运行轨迹；
8. 计算能量与频率；
9. 生成结构化结果。

当两个势阱中心或束腰不同，比较图和共享 CSV 的公共网格必须覆盖两个势阱各自的 `center + [min,max] * waist` 范围。

### `visualization.py`

至少实现：

```python
plot_potential_comparison(...)
plot_force_comparison(...)
plot_single_trajectory(...)
plot_energy_evolution(...)
plot_frequency_comparison(...)
```

要求英文坐标和标题以避免中文字体缺失；至少 200 dpi；使用非交互后端；不调用 `plt.show()`；保存后关闭 figure。

### `io_utils.py`

- 安全创建输出目录；
- 明确拒绝任何解析后位于 `sources/` 内的输出路径；
- 保存 CSV、JSON、YAML 和 Markdown；
- 将 NumPy 类型转换为可序列化类型；
- JSON 中不得出现 NaN 或 Infinity；
- 保存程序化图像检查结果。

## 六、能量误差定义

不要使用除以带符号 (E_0) 的 `relative_energy_change`。由于束缚态 (E_0<0)，该名称和符号容易产生误解；当 (E_0\to0) 时也会病态。

每个时刻至少保存：

\[
\delta E(t)=E(t)-E_0,
\qquad
\epsilon_E(t)=\frac{E(t)-E_0}{U_0},
\qquad
|\epsilon_E(t)|=\frac{|E(t)-E_0|}{U_0}.
\]

汇总指标至少包括：

- `max_abs_energy_error_over_depth`；
- `rms_energy_error_over_depth`；
- `peak_to_peak_energy_error_over_depth`；
- 可选的能量误差线性拟合斜率，但不得把有界振荡误称为 secular drift。

能量图应优先画 (E-E_0) 或 ((E-E_0)/U_0)，避免动能、势能和负的总能量共用尺度时掩盖误差。

## 七、输出规则

`--trap` 支持 `slm`、`aod` 和 `both`。输出必须与选择一致：

- `both`：生成双势阱比较图、共享网格和两套轨迹；
- `slm`：只生成 SLM 数据与图，不得伪造 AOD 输出；
- `aod`：只生成 AOD 数据与图，不得伪造 SLM 输出。

`--trap both` 默认至少生成：

```text
outputs/level0_static_trap/demo/
├── config_used.yaml
├── static_grid.csv
├── slm_trajectory.csv
├── aod_trajectory.csv
├── metrics.json
├── image_validation.json
├── summary.md
├── potential_comparison.png
├── force_comparison.png
├── slm_trajectory.png
├── aod_trajectory.png
├── slm_energy.png
├── aod_energy.png
└── frequency_comparison.png
```

`static_grid.csv` 在 `both` 模式至少包含：

```text
x_m
x_um
slm_potential_J
slm_potential_uK
slm_force_N
aod_potential_J
aod_potential_uK
aod_force_N
```

单势阱模式只保留公共坐标列和相应势阱列。

轨迹 CSV 至少包含：

```text
time_s
time_us
position_m
displacement_from_center_m
position_um
velocity_m_per_s
acceleration_m_per_s2
kinetic_energy_J
potential_energy_J
total_energy_J
energy_error_J
energy_error_over_depth
abs_energy_error_over_depth
```

`metrics.json` 对每个选中势阱至少保存：

- 输入参数、单位、来源和假设标签；
- 平衡位置、中心偏差和力残差；
- 解析曲率、两个差分步长的数值曲率和收敛差异；
- 解析小振幅频率、数值曲率频率、非线性积分频率、轨迹估计频率；
- FFT 分辨率和有效完整周期数；
- 数值离散化误差与有限振幅频移，二者不得混为一个“相对误差”；
- 初始总能量、束缚状态；
- 最大、RMS 和峰峰值归一化能量误差；
- 最大位移和最大速度；
- 关键验收项的通过/失败状态。

`summary.md` 用简洁中文说明参数、假设、主要结果、误差分类、验证状态及未包含的后续物理。

## 八、图片验收：不得假报视觉检查

代码 agent 是否能看图取决于运行环境。必须区分程序化检查和真正视觉检查。

### 必做的程序化检查

为每张 PNG 自动检查并记录到 `image_validation.json`：

- 文件存在且文件大小大于合理下限；
- 能被 Pillow 或 Matplotlib 成功重新打开；
- 像素宽高符合预期且不为零；
- 图像不是全白、全黑或近乎单色，像素方差大于合理阈值；
- 数组数据有限，绘图前后没有 NaN/Infinity；
- x/y 数据范围非零；
- 坐标范围覆盖实际数据；
- 图例标签数与所画数据系列相符；
- figure 保存后已关闭。

这些检查只能证明文件有效且不明显为空，不能证明文字没有重叠、曲线易读或布局美观。

### 条件式视觉检查

- 若运行环境明确提供图像查看/渲染工具，逐张打开最终 PNG 并记录 `visual_review_performed: true` 以及实际发现；
- 若没有此能力，记录 `visual_review_performed: false`，在 `summary.md` 和最终回复中明确写“已完成程序化图像检查，未进行人工或模型视觉审阅”；
- 未真正打开图片时，严禁声称“已确认坐标轴合理、曲线清晰可读”或类似结论。

## 九、测试要求

使用 pytest，至少覆盖：

1. 中心势能为 (-U_0)；
2. 中心受力为零；
3. 中心两侧受力指向中心；
4. `gaussian_force` 与势能有限差分梯度一致；
5. 非零中心时数值平衡位置接近配置中心，且力残差足够小；
6. (h) 和 (h/2) 的数值曲率收敛到解析曲率；
7. 数值曲率频率接近解析小振幅频率；
8. 使用 `initial_displacement_waist: 0.01` 时，轨迹频率接近小振幅频率；
9. 使用 0.10 或更大有限振幅时，不把物理非简谐频移判为积分器失败；最好与非线性周期积分结果比较；
10. Velocity-Verlet 的位置或相位误差随步长呈二阶收敛；
11. 固定物理时长并减半步长时，比较 `rms_energy_error_over_depth` 或 `peak_to_peak_energy_error_over_depth`，检查总体二阶缩放或两者均低于明确阈值；不得要求每次采样得到的“最大漂移”严格单调；
12. 配置拒绝非正、非有限、非束缚或时间步长过大的输入；
13. `--trap slm`、`aod` 和 `both` 的 CLI smoke test 均成功，并只生成相应文件；
14. `metrics.json` 可被严格 JSON 解析且不含 NaN/Infinity；
15. 所有 PNG 通过程序化图像检查。

测试容差必须注明针对哪一类误差：

- 解析与数值曲率频率：严格数值容差；
- 小振幅测试轨迹与解析频率：严格但考虑采样估计误差；
- 默认 0.10w 轨迹与小振幅频率：允许约 1% 量级、最大不超过配置的 3%，并将其报告为非简谐频移；
- 默认轨迹与非线性积分基准：使用较严格的数值容差。

## 十、执行与最终验收

按顺序完成：

1. 检查工作区、Python 和已有依赖；
2. 条件式检查 `sources/` 是否存在；
3. 创建完整项目和 `pyproject.toml`；
4. 安装 editable package；
5. 实现配置、模型、分析、积分器、输出和绘图；
6. 编写并运行测试；
7. 修复所有真实失败，不得通过删除测试或无根据放宽容差掩盖问题；
8. 从项目目录运行完整 `--trap both` demo；
9. 严格解析并检查 CSV、JSON、YAML 和 Markdown；
10. 完成所有 PNG 的程序化检查；
11. 仅在有图像查看工具时执行并声称视觉检查；
12. 检查输出未写入 `sources/`，且未修改已有参考材料；
13. 最终回复给出实际命令退出状态、测试数量、主要数值结果、输出绝对路径和未覆盖物理过程。

最终回复不得只说“已完成”。必须明确列出：

- 创建的主要模块和打包文件；
- editable 安装和 CLI 是否真实成功；
- 测试通过/失败/跳过数量；
- SLM/AOD 的解析、曲率、非线性基准和轨迹频率；
- 数值离散化误差与非简谐频移；
- 能量误差幅度，而不是笼统的“无漂移”；
- 程序化图像检查结果；
- 是否真正执行了视觉检查；
- 输出目录及主要文件；
- Level 0 尚未包含的论文关键动力学。
