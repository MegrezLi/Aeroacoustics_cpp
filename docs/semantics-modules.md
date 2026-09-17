# S4–S5：声学量与模块接口

本轮将频带、单位和计权状态写入声学接口，并把广义 α 积分器从三叶片结构后端中拆出。默认算例仍使用三叶片、每片三个模态；没有增加塔架、平台或控制器物理模型。

## 声学量

公开接口在 [acoustic_quantities.hpp](../include/acoustic_quantities.hpp)，实现位于 [quantities.cpp](../src/acoustics/quantities.cpp)。

| 类型 | 每个频带中的量 | 单位 |
| --- | --- | --- |
| `PressurePsd` | 频带内平均声压功率谱密度 | Pa²/Hz |
| `BandMeanSquarePressure` | 频带均方声压 | Pa² |
| `BandSoundPressureLevel` | 频带声压级 | dB，参考 20 μPa |
| `BandSoundPower` | 频带声功率 | W |
| `BandSoundPowerLevel` | 频带声功率级 | dB，参考 1 pW |

这些类型不能隐式互换；值、频带与 `Weighting` 一起保存，只读访问。线性量要求有限且非负，声级允许 −∞ 表示静音，拒绝 NaN 与 +∞。声功率类型只是数据和换算接口，本轮没有实现由接收点声压推算源声功率的模型。

`FrequencyBands` 显式保存中心、下限和上限，单位均为 Hz；构造时检查有限、正值、顺序及不重叠。`integrate()` 用每带平均 PSD 乘以 Hz 带宽；它不对离散频率采样点自动插值。`merge_bands()` 只合并完整、相邻的源频带，拒绝拆分、内部缺口或遗漏。总量函数仅接收线性均方声压或声功率。

```cpp
using namespace aeroacoustics;
FrequencyBands bands({{100., 80., 120.}, {150., 120., 180.}});
PressurePsd density(bands, {2e-6, 2e-6});
auto pressure = integrate(density);              // 8e-5、12e-5 Pa²
double total = total_mean_square_pressure(pressure); // 2e-4 Pa²
auto weighted = apply_a_weighting(pressure);
auto levels = pressure_levels(weighted);         // 计权状态仍为 A
```

A 计权沿用当前 OpenFAST 参考公式，按频带中心施加能量修正；宽频带下属于中心频率近似。再次给已计权结果调用 `apply_a_weighting()` 会报错。合并已计权频带时只加能量，不重新施加目标中心的计权。

## 参考频带与 TNO

默认 34 个名义中心频率仍为 10–20000 Hz。为避免四舍五入的中心导致相邻频带轻微重叠，元数据使用二进制三分之一倍频程边界：以 1000 Hz 为索引 0，第 k 带边界为 `1000 × 2^((k ± 0.5)/3)`。原名义中心保留用于公式计算及通道标签。此定义是项目的显式频带约定，不代表新增了标准认证。

`band_section_spectrum()` 返回带单位和频带的七类声源声压级；它和整机 `AcousticAggregator` 只接受这 34 个名义中心的递增子集。选取子集后，总量只代表所选频带。原 `section_spectrum()` 和标量内核仍可计算任意递增频率采样点，但不能把密集采样点当作独立频带求和，也不会因此变成窄带噪声模型。

TNO 原来的 `2ω(√r − 1/√r)` 换算保持原运算顺序，集中为 `reference_tno_bandwidth()`，返回单独的 `ReferenceTnoBandwidth` 类型。它是原模型的带宽乘数，不能当作 Hz 带宽传给 SI PSD 积分。BPM/TNO 已输出的频带声级也不能再乘一次带宽。

为兼容现有程序，`AcousticResult::power` 仍保留字段名，实际量明确为 `〈p²〉/(20 μPa)²`，不是 W。结果携带共享只读 `AcousticMetadata`；文件接收器检查频带及计权匹配。`run.json` 增加声学量、参考声压、边界和 TNO 约定。A 计权文件标题使用 dBA，未计权使用 dB；原零能量占位及 `.mask` 规则不变。

## 积分器与结构后端

| 文件 | 职责 |
| --- | --- |
| [integrator.hpp](../include/turbine/integrator.hpp)、[integrator.cpp](../src/turbine/coupling/integrator.cpp) | 自由度布局、广义 α、Newton、数值工作区及诊断 |
| [modules.hpp](../include/turbine/modules.hpp)、[modules.cpp](../src/turbine/coupling/modules.cpp) | 模块能力、支持组合校验与布局工厂 |
| [structural_adapter.cpp](../src/turbine/coupling/structural_adapter.cpp) | 当前叶片结构向通用加速度接口的适配 |
| [solver.cpp](../src/turbine/coupling/solver.cpp) | UA、BEM、结构的推进顺序和完整求解器状态 |

`DofLayout` 保存自由度名称、位移单位、加速度尺度和耦合块。一个块可以包含任意数量的相互耦合自由度；不同块只能在本次校正的外部输入冻结后相互独立。布局检查块连续覆盖全部自由度且名称唯一。当前固定塔架后端声明三个独立块，各含三个叶片模态。

`GeneralizedAlpha` 的状态、历史加速度和 Jacobian 工作区都按布局构造，步进时复用内存。3×3 系统保留原 `solve3()` 数值路径，其余尺寸使用带主元选择的稠密消元。它适合目前的小型模态系统；大规模有限元系统仍需要稀疏线性求解后端。

每步调用顺序：

```text
Solver::advance()
  GeneralizedAlpha::predict()
  Rotor::advance_airfoils()          推进一次 UA/BEM 历史
  Rotor::evaluate_into()             预测状态下的气动力
  FixedBaseAcceleration             引用冻结载荷，准备叶片基底
  GeneralizedAlpha::correct()
    AccelerationOperator::evaluate() 残差及 Jacobian 扰动，可重复调用
      Rotor::structural_acceleration()
    Newton 校正并提交状态
```

`AccelerationOperator::evaluate()` 只接收当前块的试算状态，必须覆盖全部输出，不能推进模块历史。将来加入共享的塔架或平台自由度时，必须把相互影响的自由度放入同一耦合块，并实现对应质量、载荷、运动传递及加速度方程。仅增加 `DofDescriptor` 或放宽输入开关不能得到新物理模型。

当前模块能力表说明了稳态风、固定塔架叶片、稳态 BEM、Minnema/Pierce UA、广义 α 和声学模块的输入、输出及历史状态。`configure_modules()` 集中校验支持的组合，`structural_dof_layout()` 校验所选模块并生成布局；它们不是任意模型的插件加载器。Case 与 TurbineModel 构建均执行校验，尚未实现的配置会在求解前被拒绝。

## 状态与输出兼容

`Solver::state()` 继续提供当前叶片后端的数组视图；`generalized_state()` 与 `layout()` 提供通用状态和名称。`StepView::generalized` 在标准仿真路径中有效，借用期限同其他帧数据；原三字段 `StepView` 初始化仍可用于当前固定后端。`dynamics.csv` 根据布局生成列，默认列名、顺序及数值不变。`run.json` 新增模块配置、自由度单位/尺度及耦合块。

`SolverOptions` 的同名成员访问不变，通用 Newton 参数移至基类 `NewtonOptions`。依赖聚合位置初始化的外部代码应改用具名成员赋值。自行构造 `AcousticResult` 并交给 `FileOutput` 时须提供匹配的声学元数据。

复制 `GeneralizedAlpha` 保存它自身的状态和算法历史，不包含外部物理模块。整机继续使用 `Solver`/`Simulation` 的完整检查点，涵盖 UA、BEM、TI 与采样进度；失败后须恢复或重置，禁止直接推进。

## 验证

[semantics_modules_probe.cpp](../tests/semantics_modules_probe.cpp) 验证已知 PSD 积分、三分之一倍频程合并、声压/声功率换算、A 计权及非法频带；另用含非对角耦合的 1、2、3、4、6 自由度系统核对解析解，并检查复制恢复、重置及模块组合拒绝。

[check_semantics.py](../tests/check_semantics.py) 检查实际 CLI 的频带、单位、dB/dBA 标题、自由度元数据及未支持模块报错。Python 只组织运行和检查输出。

```sh
./build/semantics_modules_probe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst
python tests/check_semantics.py build/aeroacoustics_turbine build/semantics-check
```

本轮同时保留声学 Fortran 对照、8/9 m/s 整机对照、结构零分配检查、检查点及并行回归。验证记录见 [validation-semantics.json](validation-semantics.json)。
