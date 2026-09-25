# 开发接口与模块

[源码索引](#source) · [参考工况与限制](#scope) · [仿真与状态](#simulation) · [声学量与积分器](#quantities-modules) · [工具链](#toolchain)

从输入到函数调用的完整顺序见 [工作流](../aeroacoustics工作流.md)，闭环控制与传播模型见 [工程模型](engineering.md)，测试命令集中在 [验证](validation.md)。

<a id="source"></a>

## 源码索引

源码按物理模型和程序职责组织。公开接口位于仓库根目录的 `include/`，命令行程序入口位于 `src/apps/`。

| 目录 | 文件 | 职责 |
| --- | --- | --- |
| `acoustics/` | [boundary_layer.cpp](../src/acoustics/boundary_layer.cpp) | BPM 经验边界层、表格边界层参数转换 |
| | [bpm_trailing_edge.cpp](../src/acoustics/bpm_trailing_edge.cpp) | 压力面、吸力面及分离尾缘噪声 |
| | [bpm_other_sources.cpp](../src/acoustics/bpm_other_sources.cpp) | 层流、尾缘钝度和叶尖噪声 |
| | [bpm_spectral_shapes.cpp](../src/acoustics/bpm_spectral_shapes.cpp)、[directivity.cpp](../src/acoustics/directivity.cpp) | BPM 分段谱形函数与声源指向性 |
| | [inflow_noise.cpp](../src/acoustics/inflow_noise.cpp) | Lowson 入流噪声和 Simplified Guidati 修正 |
| | [tno.cpp](../src/acoustics/tno.cpp) | TNO 尾缘噪声与模型积分装配 |
| | [spectrum.cpp](../src/acoustics/spectrum.cpp) | 参数检查、声源准备与观察点频谱装配、A 计权和声能叠加工具 |
| | [quantities.cpp](../src/acoustics/quantities.cpp) | 频带边界、PSD 积分、声压/声功率换算、计权与频带合并 |
| | [workspace.cpp](../src/acoustics/workspace.cpp)、[source_models.hpp](../src/acoustics/source_models.hpp) | 跨观察点共享声源计算，复用频谱与 TNO 积分缓冲区；内部模型数据类型 |
| | [kernel_compat.cpp](../src/acoustics/kernel_compat.cpp) | 原标量函数接口的适配，用于直接调用及 Fortran 公式对照 |
| | [geometry.cpp](../src/acoustics/geometry.cpp) | 观察点与叶片前后缘坐标变换 |
| | [driver.cpp](../src/acoustics/driver.cpp) | 声学时间状态、边界层表、湍流强度和声学输入读取 |
| `numerics/` | [quadrature.cpp](../src/numerics/quadrature.cpp) | 61 点 Gauss–Kronrod 积分及误差估计 |
| | [backend.cpp](../src/numerics/backend.cpp) | 普通/MKL 后端标识与向量点积 |
| `propagation/` | [outdoor.cpp](../src/propagation/outdoor.cpp) | 空气吸收、镜像地面反射、主导薄屏障衍射与频带干涉 |
| `turbine/inflow/` | [grid_wind.cpp](../src/turbine/inflow/grid_wind.cpp) | 时间/x/y/z 网格的三分量速度插值，拒绝越界 |
| `turbine/control/` | [controller.cpp](../src/turbine/control/controller.cpp) | 两质量传动链、发电机转矩、PI 变桨与偏航伺服 |
| `turbine/analysis/` | [statistics.cpp](../src/turbine/analysis/statistics.cpp) | 接收时间序列、声能积分、时间百分位、调制分析和声功率几何换算 |
| | [metrics.cpp](../src/turbine/analysis/metrics.cpp) | 节点历史、音调输运、受声点统计、风速分箱和报告输出 |
| | [validation.cpp](../src/turbine/analysis/validation.cpp) | 配对误差、有限差分敏感性、相关不确定性与等权样本统计 |
| `turbine/aerodynamics/` | [bem.cpp](../src/turbine/aerodynamics/bem.cpp)、[unsteady.cpp](../src/turbine/aerodynamics/unsteady.cpp) | BEM 诱导求解、偏斜修正和非定常翼型状态 |
| `turbine/structure/` | [blade_dynamics.cpp](../src/turbine/structure/blade_dynamics.cpp) | 叶片模态、运动学、广义载荷与加速度 |
| `turbine/coupling/` | [mesh.cpp](../src/turbine/coupling/mesh.cpp) | 气动和结构网格的运动、载荷传递 |
| | [rotor.cpp](../src/turbine/coupling/rotor.cpp) | 叶素状态、气动力和转子装配 |
| | [solver.cpp](../src/turbine/coupling/solver.cpp) | UA、BEM、结构的推进顺序和完整求解状态 |
| | [integrator.cpp](../src/turbine/coupling/integrator.cpp) | 按自由度布局执行广义 α 与 Newton 迭代 |
| | [structural_adapter.cpp](../src/turbine/coupling/structural_adapter.cpp) | 当前三叶片结构后端的加速度接口 |
| | [modules.cpp](../src/turbine/coupling/modules.cpp) | 模块能力表、组合校验和结构布局工厂 |
| `turbine/io/` | [case_input.cpp](../src/turbine/io/case_input.cpp) | 输入解析、配置检查、翼型表及稳态风场 |
| | [engineering.cpp](../src/turbine/io/engineering.cpp) | 独立传播配置与屏障表读取 |
| | [surface.cpp](../src/turbine/io/surface.cpp)、[metrics_input.cpp](../src/turbine/io/metrics_input.cpp) | 表面数据及工程统计配置读取、来源和范围检查 |
| | [validation_input.cpp](../src/turbine/io/validation_input.cpp) | 独立数据 CSV、分组/上下文检查、接收点导入和验证报告 |
| `interfaces/` | [c_api.cpp](../src/interfaces/c_api.cpp) | 对外 C ABI、错误与调用状态 |
| `turbine/simulation/` | [simulation.cpp](../src/turbine/simulation/simulation.cpp) | 完整仿真调度、重置和检查点 |
| | [acoustic_adapter.cpp](../src/turbine/simulation/acoustic_adapter.cpp) | 声学配置装配、气动状态到声学节点的转换 |
| | [results.cpp](../src/turbine/simulation/results.cpp) | 按输出需求分配通道、分块聚合声能 |
| | [batch.cpp](../src/turbine/simulation/batch.cpp) | 独立工况线程调度、输出路径冲突检查和逐项错误 |
| | [file_output.cpp](../src/turbine/simulation/file_output.cpp) | 标准文件输出及失败状态 |
| `apps/` | [section_main.cpp](../src/apps/section_main.cpp) | 单截面频谱示例入口 |
| | [turbine_main.cpp](../src/apps/turbine_main.cpp) | 整机命令行参数解析及运行入口 |
| | [batch_main.cpp](../src/apps/batch_main.cpp) | 批量工况命令行入口 |
| | [validation_main.cpp](../src/apps/validation_main.cpp) | 独立验证与不确定性命令行入口 |

`apps/turbine_main.cpp` 调用 `run_case()`；`Simulation` 通过适配器、声学驱动及聚合器处理各时间步，结果交给 `ResultSink`。整机求解器负责协调气动、结构与网格传递；声学模型调用数值积分和后端计算。

声学主路径为 `AcousticDriver::step_blocks()` → `AcousticWorkspace::evaluate_blocks()` → `AcousticAggregator::append()`。每个采样时刻先对各节点调用 `prepare_section()`，再按观察点调用 `emit_section()`。各模型的 `prepare_*()` 计算与观察点无关的谱形，`emit_*()` 施加距离和指向性。源码调用链见 [工作流](../aeroacoustics工作流.md)，性能和数值检查见 [P4–P7 说明](performance.md#p4-p7)。完整快照接口仍保留。

可靠性辅助接口位于 `include/`：[acoustic_levels.hpp](../include/acoustic_levels.hpp) 区分静音与非法声级；[checked_output.hpp](../include/checked_output.hpp) 检查输出生命周期；[lookup_diagnostics.hpp](../include/lookup_diagnostics.hpp) 提供显式、线程局部的查表诊断会话。掩码、报告和严格模式见 [R1–R3 说明](validation.md#reliability)。

### 构建文件

- 根目录 `CMakeLists.txt`：构建开关、编译标准、MKL 查找和输出位置。
- [CMakeLists.txt](../src/CMakeLists.txt)：声学静态库及动态库的共用源文件清单。
- [turbine/CMakeLists.txt](../src/turbine/CMakeLists.txt)：整机模块库。
- [apps/CMakeLists.txt](../src/apps/CMakeLists.txt)：截面、整机、批量和独立验证四个可执行程序。
- 根目录 `tests/CMakeLists.txt`：数值对照探针。

公开头文件路径、构建目标名称及可执行文件位置保持兼容。原有 `section_spectrum()`、`snapshot_spectrum()` 和返回独立快照的 `AcousticDriver::step()` 仍可使用；连续计算可使用复用缓冲区的接口。构建及运行命令见根目录 README。

<a id="scope"></a>

## 参考工况与支持范围

本文描述官方固定工况的参考路径。新增非稳态风、闭环传动链/执行器和室外传播通过独立配置启用，见 [E1–E2 工程说明](engineering.md)；其简化物理模型与原 Fortran 对照范围分开记录。

`aeroacoustics_turbine` 的运行入口是 `src/apps/turbine_main.cpp`。CMake 仅启用 CXX，默认构建不链接 OpenFAST、Fortran、BLAS 或 MKL。输入来自官方算例文件，风机状态从初始条件逐步计算。

表内 `.cpp` 路径相对于 `src/`。

| C++ 文件 | 对应模块与计算内容 |
| --- | --- |
| `turbine/io/case_input.cpp` | 输入字段、翼型与坐标文件、AirfoilInfo 线性插值、InflowWind 恒定风 |
| `turbine/structure/blade_dynamics.cpp` | ElastoDyn 叶片三模态方程、运动学、广义质量与载荷 |
| `turbine/aerodynamics/bem.cpp` | BEMTUncoupled 诱导公式、工作区间及入流角求根、偏斜修正 |
| `turbine/aerodynamics/unsteady.cpp` | UA_Mod=3 的附着流、分离、涡升力、滤波及离散状态 |
| `turbine/coupling/mesh.cpp`、`include/turbine/math.hpp` | NWTC 运动与载荷映射、旋转插值及坐标变换 |
| `turbine/coupling/rotor.cpp` | 转子盘与叶素坐标、风速、气动力、UA/BEM 状态传递 |
| `turbine/coupling/solver.cpp` | 气动／结构调用顺序和完整状态管理 |
| `turbine/coupling/integrator.cpp`、`structural_adapter.cpp` | 通用广义 α 积分、Newton 及叶片结构适配 |
| `turbine/coupling/modules.cpp` | 支持的模块组合、能力声明与自由度布局 |
| `turbine/simulation/` | 仿真调度、声学节点适配、声能聚合及结果输出 |
| `apps/turbine_main.cpp` | 命令行参数解析和库入口调用 |

### 官方配置

3 个叶片，叶片柔性长度 63 m，轮毂半径 2 m。每片有 17 个结构积分节点，加上根尖运动节点；气动网格为 30 个预弯节点。每片启用两个挥舞自由度和一个摆振自由度，共 9 个自由度。

结构方程包括模态形状及结构扭角、轴向缩短、质量矩阵、刚度、阻尼、重力、离心力和科里奥利项。转子以 10.04 rpm 匀速旋转，固定桨距 1.17°、预锥 −2.5°、轴倾 −5°。塔架和其余自由度锁定，与官方输入一致。

气动使用 `BEM_Mod=1`、`Wake_Mod=1`、`DBEMT_Mod=0`、`UA_Mod=3`、`FLookup=True`、`AFTabMod=1`。上游固定版本的 AirfoilInfo 默认插值阶数为 **1**，输入注释中的“默认 3”不适用于该版本。

### 每步求解顺序

1. 依据广义 α 格式预测叶片模态位移、速度及方位角。
2. 将预测运动映射至气动网格，用上一时刻输入推进非定常翼型状态。
3. 计算当前 BEM、气动力与声学输入。BEM 保留上一步的入流角，按原工作区间缩小括区，采用 Brent–Dekker 求根。
4. 保持该步气动载荷，在结构迭代中更新载荷作用点映射，求出结构加速度并校正位移、速度。
5. 按 `DT_AA` 输出声学结果；结构状态按每个 `DT` 输出。

这种调用顺序对应原 FAST_Solver 的 AeroDyn Option 2 和紧耦合 ElastoDyn。初始物理及算法加速度均设为零，与原求解器 Step0 后的状态一致。

默认参考模式下，C++ 将固定塔架条件下的结构方程分解为三个 3×3 系统，每次 Newton 迭代重新计算数值 Jacobian，修正量阈值为 1e−9，最多 12 次。原 FAST_Solver 使用包含模块输入的整体 Jacobian、缓存更新及不同的停止准则。因此 `ConvTol`、`MaxConvIter`、`DT_UJac`、`UJacSclFact` 不控制本实现的内部 Newton 迭代。`RhoInf`、`DT` 参与广义 α 系数计算；要求 `ModCoupling=3`、`NumCrctn=0`。

可通过 `SolverOptions` 修改内部迭代参数，或启用 `scaled` 模式的尺度化差分、残差判停和同一步内 Jacobian 复用；详见 [P7 求解说明](performance.md#p4-p7)。

每次 Jacobian 扰动只改变一片叶片，使用 `Rotor::structural_acceleration()` 共享该状态的节点运动，更新载荷映射并计算加速度。每轮 Newton 迭代的单叶片映射次数从 30 次降为 12 次，积分公式、扰动步长和停止准则不变。偏斜模型开关及系数在 `Rotor` 初始化时缓存。

声学驱动通过 `step_blocks()` 回调逐块交付观察点频谱，由 `AcousticWorkspace` 持有当前块，`AcousticAggregator::append()` 累计请求的输出。完整 `step_view()` 接口仍保留。工作区先准备各节点的边界层和声源谱形，再分别计算观察点的距离与指向性；主程序汇总输出后才进入下一步。频谱、TNO 积分和输出数组反复使用，湍流强度仍按原来的每步更新顺序推进。详见 [工作流](../aeroacoustics工作流.md) 和 [优化记录](performance.md#initial)。

BEM 使用 `IndToler` 和 `MaxIter`，默认双精度残差阈值 5e−10，并使用 1e−6 rad 的括区停止阈值。其插值和停止路径与原 Fortran Brent 例程不同，未复制原 `mod_root1dim.f90`。这些数值实现差异在整机误差报告中保留。

### 输入限制

实现覆盖上面的官方配置及已验证的 9 m/s 风速变体。支持固定转速、桨距、稳态风及表格参数的输入读取，但两项回归通过不代表任意参数组合都已验证。

主要限制由 `Case::validate_scope` 及各模块构造函数检查：

模块与自由度组合校验集中在 `configure_modules()`，Case 和冻结后的 TurbineModel 均执行检查。数值积分器已支持不同尺寸的耦合块，但当前整机物理后端仍仅实现下列范围，详见 [S4–S5](development.md#quantities-modules)。

- 只支持单转子、三叶片、相同结构与气动叶片文件；叶片三个模态必须全部启用。
- 仅支持 WindType=1、单张翼型表和线性翼型插值。
- 原 OpenFAST 输入中的塔架、平台、传动链、变桨、偏航动力学及控制器开关保持关闭；独立工程控制器通过 `--controller` 配置。固定平台偏移、初始叶片挠度、叶尖附加质量要求为零。
- 不支持动态入流、自由涡尾流、其他 UA 模型、塔影、塔架／机舱／尾翼气动力或扇区平均。
- 气动步长与主步长相同，声学步长须为其整数倍；整机驱动支持 `TICalcMeth=1`。
- 不提供 OpenFAST 的线性化、稳态求解、重启文件或 `.outb` 输出；C++ 库提供独立的进程内完整检查点，见 [仿真接口](development.md#simulation)。`OutList`、`DT_Out` 等原动力学输出设置不用于选择 C++ 的 `dynamics.csv` 通道。

声学库还含 TNO、层流与叶尖模型以及其他辅助算法，默认整机案例的开关并不覆盖所有声学分支。源码级 Fortran 对照单独验证这些公式。

<a id="simulation"></a>

## 仿真接口、声源配置与状态

整机调度、声源配置与状态管理分别提供公开接口。

### S1：仿真与输出分开

`main()` 负责解析命令行参数、调用 `run_case()`，并打印完成信息或错误。库接口不依赖命令行局部变量，也不在求解过程中直接打印日志。

| 层 | 接口 | 职责 |
| --- | --- | --- |
| 输入与模型 | `Case`、`TurbineModel` | 读取输入，建立不可变的结构和气动模型数据 |
| 声学装配 | `AcousticConfiguration` | 集中取得介质参数、观察点、采样间隔和 TI 配置 |
| 数据适配 | `AcousticInputAdapter` | 将 `RotorOutput` 转为声学节点，按采样计划插值边界层 |
| 时间调度 | `Simulation::next()` | 协调结构求解、声学采样、TI 更新和结果聚合 |
| 结果汇总 | `AcousticAggregator`、`OutputLayout` | 按观察点、频带、机制和节点累加声能，建立通道标签 |
| 输出 | `ResultSink`、`FileOutput` | 接收帧数据；默认实现写入原 `.out`、`.mask`、CSV 和元数据 |

这些实现位于 [src/turbine/simulation](../src/turbine/simulation)。`AcousticResult::power` 保存相对均方声压；`FileOutput` 写文件时才转为 dB，零声能占位和异常声级处理沿用 R1–R3。

现有命令行保持不变。其他 C++ 程序可以链接 `turbine_modules`，直接使用：

```cpp
#include "turbine/simulation.hpp"

turbine::RunOptions options;
options.duration = 2.0;
auto summary = turbine::run_case("case.fst", "results", options);
```

不写文件时，逐帧读取结果：

```cpp
turbine::TurbineModel model{turbine::Case("case.fst")};
turbine::Simulation simulation(model, options);
while (auto frame = simulation.next()) {
    // frame->state 是只读结构状态。
    if (frame->acoustics) {
        const auto &observer_power = frame->acoustics->power[0];
        // 在这里交给优化目标函数、统计程序或其他结果接收器。
    }
}
auto summary = simulation.summary();
```

也可以实现 `ResultSink::begin/write/finish`，再调用 `run(simulation, sink)`。该辅助函数要求全新或已重置的仿真，防止把半程结果当作完整文件写出；分段运行使用 `next()`。

`StepView`、节点适配结果和聚合结果均为借用视图，下一次推进、恢复、重置、赋值或对象销毁后失效。需要长期保存时，调用者应复制所需数据。默认时间循环仍复用工作区。

### S2：机制与模型分别描述

[mechanisms.hpp](../include/mechanisms.hpp) 集中定义 `Mechanism` 枚举、通道名称、单位、顺序及是否计入总声能。机制数组长度、输出标签、聚合循环和诊断名称都从这份描述表取得。

| 输入选择 | 内部模型选择 | 通道处理 |
| --- | --- | --- |
| `TBLTEMod=0` | `TrailingEdgeModel::off` | 三个尾缘通道关闭 |
| `TBLTEMod=1` | `TrailingEdgeModel::bpm` | 压力面、吸力面、分离均使用 BPM |
| `TBLTEMod=2` | `TrailingEdgeModel::tno_with_bpm_separation` | TNO 替代压力面和吸力面，分离仍使用 BPM |
| `TIMod=1` | `InflowModel::lowson` | Lowson 入流噪声 |
| `TIMod=2` | `InflowModel::lowson_guidati` | Lowson 加 Guidati 修正，仍为同一入流通道 |

`SourceSelection` 在声学工作区初始化时解析，准备和发射阶段使用同一份选择。TNO 和 BPM 的完整尾缘结果不会被叠加计算。机制仍按 `LBL / TBL_pressure / TBL_suction / TBL_separation / bluntness / tip / inflow` 输出。

新增模型通常只需增加模型选择和对应准备、发射函数；新增物理机制则还需扩展枚举、描述表及通道兼容测试。现有 C ABI 的操作编号和缓冲区布局没有改动；它们是截面公式接口，不等于机制枚举。

### S3：模型所有权与完整状态

`TurbineModel` 深复制传入的 `Case`，向外只提供 const 访问。`Rotor` 共享这份模型；其非定常翼型对象通过共享所有权引用其中的翼型，不依赖调用者容器的地址。单独构造 `UnsteadyAirfoil(const Airfoil&, ...)` 时，该对象复制并持有翼型数据。

`Solver` 的结构状态、气动结果、时间、加速度和转子改为私有成员。常用读取方式改为：

| 原调用 | 当前调用 |
| --- | --- |
| `solver.state` | `solver.state()` |
| `solver.aerodynamic` | `solver.aerodynamic()` |
| `solver.time`、`solver.step_number` | `solver.time()`、`solver.step_number()` |
| `solver.rotor.structure` | `solver.rotor().structure()` |
| `driver.state.values` | `driver.turbulence_state().values` |

这些访问器返回值或只读引用，不能单独改写时间或位移，绕过相应的气动和历史状态更新。需要不同初始模型时构造新的 `Case`/`TurbineModel`。当前不支持任意初始挠度。

`Solver` 与 `Simulation` 均支持复制、移动、`checkpoint()`、`restore()` 和 `reset()`。复制后的实例共享不可变模型，演化状态及缓冲区独立。移动后的源对象仅用于销毁或重新赋值。

```cpp
auto saved = simulation.checkpoint();
auto branch = simulation;       // 独立推进的副本
branch.next();
simulation.restore(saved);      // 恢复完整状态
```

求解器检查点包括结构状态、积分器加速度、气动结果、BEM 前一步根、UA 离散历史及工作区；仿真检查点另外包含声学采样状态、TI 缓冲区、节点适配数据、聚合结果和诊断记录。只接受来自同一个 `TurbineModel` 共享实例的检查点，避免把其他模型的历史混入当前运行。

`reset()` 从已加载数据恢复初始状态，不重新打开文件。重新读取已修改的输入，应构造新的模型和仿真。`Simulation` 首次构造仍会读取声学观察点和边界层文件；构造完成后运行、复制、恢复和重置不依赖这些文件继续存在。

求解异常后，实例拒绝继续推进或建立检查点，须先恢复或重置。文件写入异常后，`FileOutput` 拒绝继续写入成功元数据，须创建新的输出器。检查点为进程内对象，不是 OpenFAST 重启文件，也不包含输出文件的回滚或追加；分段输出需要调用者管理接收器。

`run.json` 的物理时间、步数和采样数保持不变；`elapsed_seconds` 现在统计从 `Simulation` 初始化到最后一帧处理完毕的墙钟时间，不含此前的 `Case` 文件解析及最后的输出关闭。跨版本性能比较应使用外部墙钟计时。

### P4–P7 接口补充

`RunOptions` 新增 `observer_block_size`（默认 1）和 `solver`（`SolverOptions`）。`Simulation` 的声学路径改用观察点分块回调，`AcousticResult::power` 中未由 `NrOutFile` 请求的类别为空。需要完整四类结果时设置 `NrOutFile=4`；已有完整快照 API 不变。

`Solver::diagnostics()` 与 `RunSummary::structure` 提供求解统计。默认参考模式保持原数值路径，尺度化模式须明确选择。`run.json` 新增 `structural_*` 字段。`SolverOptions` 随检查点恢复，`reset()` 保留选项。

`turbine/batch.hpp` 提供独立工况的 `run_cases()`。线程、路径限制、借用块生命周期、完整参数与验证见 [P4–P7 说明](performance.md#p4-p7)。

### S4–S5 接口补充

`AcousticResult` 新增只读共享元数据；`power` 表示相对均方声压，不是瓦特。`FileOutput` 要求声学频带与计权匹配。A 计权标题现在使用 dBA。

`Solver::generalized_state()`、`layout()` 提供通用状态与自由度布局，`StepView::generalized` 为对应借用视图。`Solver::state()` 和原三字段 StepView 用法仍适用于固定后端；输出列改由布局生成。`SolverOptions` 的 Newton 参数移至 `NewtonOptions` 基类，成员名不变，建议使用具名成员赋值。

广义 α 与 Newton 独立为 `GeneralizedAlpha`，模块组合校验集中到 `configure_modules()`。整机完整检查点范围不变。频带转换、结构适配和新增物理后端的要求见 [声学量与模块接口](development.md#quantities-modules)。

### E1–E2 库配置

`SolverOptions::wind` 可注入共享只读 `WindField`，`controller` 可提供 `ControlConfig`；`RunOptions::propagation` 提供传播参数。风场实现的 const 查询必须线程安全，不得隐藏演化状态。默认均为空，保留参考工况。

`StepView::operation` 在启用控制器时借用其 `ControlState`，包含方位角、转速、发电机转矩、轴扭转、桨距与偏航等；有效期同其他帧引用。完整检查点复制所有控制状态；共享风场与传播配置保持只读。九个叶片模态的 `generalized_state()` 不包含独立传动链/执行器的一阶状态。配置例子及耦合近似见 [E1–E2](engineering.md)。

### 表面数据与工程统计接口（E3、E5、E7）

```cpp
turbine::RunOptions options;
options.surfaces = turbine::SurfaceSet::read("surfaces.dat");
options.metrics = turbine::MetricsOptions::read("metrics.dat");
auto result = turbine::run_case("case.fst", "results", options);
```

`SurfaceSet` 只读共享，包含已加载的完整数据。`Simulation` 初始化时应用极曲线覆盖并建立独立的冻结模型；叶素适配器只对有状态数据的截面使用指定边界层。原未覆盖截面的设置不变。

`EngineeringMetrics` 独立于原四类聚合结果，在声学分块回调中记录每个节点的接收历史。历史随仿真复制/检查点深复制，`reset()` 清空；失败后禁止继续追加或导出，须恢复有效副本或重新构造。`RunSummary::metrics` 在启用时提供独立只读历史快照，可调用 `write(directory)` 输出报告；目录须已存在。默认文件接收器自动完成此步，并在所有文件成功关闭后写入 `run.json`。

历史内存按观察点、节点、频带和采样数增长；`MaxValues` 限制数值载荷个数，不包括容器容量与对象开销。仿真复制、检查点及结果快照也会复制历史，不能把这个上限当作进程 RSS 上限。长时间/大地图任务应明确选择采样间隔与数据量预算。

`ArrivalSeries`、`level_statistics()`、`modulation()` 和 `apparent_sound_power()` 可单独用于分析。它们接收具名输入并检查范围，具体时间、量纲和简化假设见 [工程模型](engineering.md#receiver-metrics)。原七类 `.out` 不包含外部音调，接收统计文件单独包含这些输入声源。

<a id="quantities-modules"></a>

## 声学量、积分器与物理后端

### 声学量

公开接口在 [acoustic_quantities.hpp](../include/acoustic_quantities.hpp)，实现位于 [quantities.cpp](../src/acoustics/quantities.cpp)。

| 类型 | 每个频带中的量 | 单位 |
| --- | --- | --- |
| `PressurePsd` | 频带内平均声压功率谱密度 | Pa²/Hz |
| `BandMeanSquarePressure` | 频带均方声压 | Pa² |
| `BandSoundPressureLevel` | 频带声压级 | dB，参考 20 μPa |
| `BandSoundPower` | 频带声功率 | W |
| `BandSoundPowerLevel` | 频带声功率级 | dB，参考 1 pW |

这些类型不能隐式互换；值、频带与 `Weighting` 一起保存，只读访问。线性量要求有限且非负，声级允许 −∞ 表示静音，拒绝 NaN 与 +∞。声功率类型只是数据和换算接口，目前没有实现由接收点声压推算源声功率的模型。

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

### 参考频带与 TNO

默认 34 个名义中心频率仍为 10–20000 Hz。为避免四舍五入的中心导致相邻频带轻微重叠，元数据使用二进制三分之一倍频程边界：以 1000 Hz 为索引 0，第 k 带边界为 `1000 × 2^((k ± 0.5)/3)`。原名义中心保留用于公式计算及通道标签。此定义是项目的显式频带约定，不代表新增了标准认证。

`band_section_spectrum()` 返回带单位和频带的七类声源声压级；它和整机 `AcousticAggregator` 只接受这 34 个名义中心的递增子集。选取子集后，总量只代表所选频带。原 `section_spectrum()` 和标量内核仍可计算任意递增频率采样点，但不能把密集采样点当作独立频带求和，也不会因此变成窄带噪声模型。

TNO 原来的 `2ω(√r − 1/√r)` 换算保持原运算顺序，集中为 `reference_tno_bandwidth()`，返回单独的 `ReferenceTnoBandwidth` 类型。它是原模型的带宽乘数，不能当作 Hz 带宽传给 SI PSD 积分。BPM/TNO 已输出的频带声级也不能再乘一次带宽。

为兼容现有程序，`AcousticResult::power` 仍保留字段名，实际量明确为 `〈p²〉/(20 μPa)²`，不是 W。结果携带共享只读 `AcousticMetadata`；文件接收器检查频带及计权匹配。`run.json` 增加声学量、参考声压、边界和 TNO 约定。A 计权文件标题使用 dBA，未计权使用 dB；原零能量占位及 `.mask` 规则不变。

### 积分器与结构后端

| 文件 | 职责 |
| --- | --- |
| [integrator.hpp](../include/turbine/integrator.hpp)、[integrator.cpp](../src/turbine/coupling/integrator.cpp) | 自由度布局、广义 α、Newton、数值工作区及诊断 |
| [modules.hpp](../include/turbine/modules.hpp)、[modules.cpp](../src/turbine/coupling/modules.cpp) | 模块能力、支持组合校验与布局工厂 |
| [structural_adapter.cpp](../src/turbine/coupling/structural_adapter.cpp) | 当前叶片结构向通用加速度接口的适配 |
| [solver.cpp](../src/turbine/coupling/solver.cpp) | UA、BEM、结构的推进顺序和完整求解器状态 |

`DofLayout` 保存自由度名称、位移单位、加速度尺度和耦合块。一个块可以包含任意数量的相互耦合自由度；不同块只能在本次校正的外部输入冻结后相互独立。布局检查块连续覆盖全部自由度且名称唯一。当前固定塔架后端声明三个独立块，各含三个叶片模态。

`GeneralizedAlpha` 的状态、历史加速度和 Jacobian 工作区都按布局构造，步进时复用内存。3×3 系统保留原 `solve3()` 数值路径，其余尺寸使用带主元选择的稠密消元。它适合目前的小型模态系统；大规模有限元系统仍需要稀疏线性求解后端。

未开启工程控制器时，每步调用顺序如下；开启后的控制/传动链交换见 [工程模型](engineering.md)。

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

### 状态与输出兼容

`Solver::state()` 继续提供当前叶片后端的数组视图；`generalized_state()` 与 `layout()` 提供通用状态和名称。`StepView::generalized` 在标准仿真路径中有效，借用期限同其他帧数据；原三字段 `StepView` 初始化仍可用于当前固定后端。`dynamics.csv` 根据布局生成列，默认列名、顺序及数值不变。`run.json` 新增模块配置、自由度单位/尺度及耦合块。

`SolverOptions` 的同名成员访问不变，通用 Newton 参数移至基类 `NewtonOptions`。依赖聚合位置初始化的外部代码应改用具名成员赋值。自行构造 `AcousticResult` 并交给 `FileOutput` 时须提供匹配的声学元数据。

复制 `GeneralizedAlpha` 保存它自身的状态和算法历史，不包含外部物理模块。整机继续使用 `Solver`/`Simulation` 的完整检查点，涵盖 UA、BEM、TI 与采样进度；失败后须恢复或重置，禁止直接推进。

<a id="toolchain"></a>

## 工具链与可选 MKL

默认独立 C++ 程序需要 CMake ≥3.20 和 C++17 编译器。GFortran 只用于原源码对照；MKL 为可选声学后端。

### 本机 GCC

已安装 GCC / G++ / GFortran 16.2.0，来自 [WinLibs](https://winlibs.com/)，x86_64、POSIX、SEH、UCRT，MinGW-w64 14.0.0。

- 发布标签：`16.2.0posix-14.0.0-ucrt-r1`。
- 安装位置：`D:/Code_Configuration/gcc-16.2.0/mingw64/bin`。
- 安装包 SHA256：`c1f52294597c0b73786b2a78eb5d176d89226d2f21875eab75e783a8b1cefcc4`，已与发布文件校验值核对。
- 原 GCC 7.3 保留在 `D:/Code_Configuration/mingw64/bin`。
- 验证使用 CMake 4.1.0；原 Fortran 的 BLAS/LAPACK 使用 Intel oneMKL 2025.2。

`. ./tools/use_gcc.ps1` 设置当前终端的编译器及 PATH，不修改永久系统配置。只有执行 Fortran 对照时才需要 GFortran。

### 可选 MKL 声学后端

```powershell
. ./tools/use_gcc.ps1
$mklRoot = 'C:/Program Files (x86)/Intel/oneAPI/mkl/2025.2'
$env:PATH += ";$mklRoot/bin"
$env:AEROACOUSTICS_RUNTIME_DIRS += ";$mklRoot/bin"
$env:MKL_THREADING_LAYER = 'SEQUENTIAL'
cmake -S . -B build-mkl -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_USE_MKL=ON "-DMKL_DIR=$mklRoot/lib/cmake/mkl"
cmake --build build-mkl -j 4
./build-mkl/aeroacoustics_example.exe build-mkl/spectrum.csv
python tests/compare_fortran.py build-mkl/libaeroacoustics_shared.dll --fortran-compiler $env:FC --report build-mkl/kernels.json
```

MinGW 通过 `mkl_rt` 单一动态接口链接；运行时应能找到 MKL DLL。TNO 积分调用 `vdExp`、`cblas_dgemv` 和 `cblas_ddot`，是否加速取决于工况和规模。

### 原 Fortran 整机参考

下载固定原源码后，使用同一 PowerShell 会话：

```powershell
cmake -S _work/openfast -B _work/build-reference -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DDOUBLE_PRECISION=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 "-DBLAS_LIBRARIES=$mklRoot/lib/mkl_rt.lib" "-DLAPACK_LIBRARIES=$mklRoot/lib/mkl_rt.lib" -DBLAS_FOUND=TRUE -DLAPACK_FOUND=TRUE
cmake --build _work/build-reference --target openfast -j 4
python tools/run_full_case.py --executable _work/build-reference/glue-codes/openfast/openfast.exe --label fortran-reference
```

OpenFAST 双精度配置使用 `-fdefault-real-8 -fdefault-double-8`。原 Fortran 参考链接 MKL；独立 C++ 默认构建不链接 MKL，也不复用这个 OpenFAST 可执行文件。两者的差异保存在验证报告中。


E8 公开接口见 [validation.hpp](../include/turbine/validation.hpp)，位于 `turbine::validation` 命名空间。`residual/error_statistics`、`uncertainty_budget`、`ensemble_statistics` 均为不依赖求解器状态的数值函数；CSV、配对规则和文件输出在独立 IO 源文件中。所有函数按值返回结果，无全局可变状态。后续可从仿真结果或外部数据调用这些接口，保持数据准备、模型运行与统计分析分离。
