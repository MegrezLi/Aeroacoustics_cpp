# 源码目录

源码按物理模型和程序职责组织。公开接口位于仓库根目录的 `include/`，命令行程序入口位于 `apps/`。

| 目录 | 文件 | 职责 |
| --- | --- | --- |
| `acoustics/` | [boundary_layer.cpp](acoustics/boundary_layer.cpp) | BPM 经验边界层、表格边界层参数转换 |
| | [bpm_trailing_edge.cpp](acoustics/bpm_trailing_edge.cpp) | 压力面、吸力面及分离尾缘噪声 |
| | [bpm_other_sources.cpp](acoustics/bpm_other_sources.cpp) | 层流、尾缘钝度和叶尖噪声 |
| | [bpm_spectral_shapes.cpp](acoustics/bpm_spectral_shapes.cpp)、[directivity.cpp](acoustics/directivity.cpp) | BPM 分段谱形函数与声源指向性 |
| | [inflow_noise.cpp](acoustics/inflow_noise.cpp) | Lowson 入流噪声和 Simplified Guidati 修正 |
| | [tno.cpp](acoustics/tno.cpp) | TNO 尾缘噪声与模型积分装配 |
| | [spectrum.cpp](acoustics/spectrum.cpp) | 参数检查、声源准备与观察点频谱装配、A 计权和声能叠加工具 |
| | [workspace.cpp](acoustics/workspace.cpp)、[source_models.hpp](acoustics/source_models.hpp) | 跨观察点共享声源计算，复用频谱与 TNO 积分缓冲区；内部模型数据类型 |
| | [kernel_compat.cpp](acoustics/kernel_compat.cpp) | 原标量函数接口的适配，用于直接调用及 Fortran 公式对照 |
| | [geometry.cpp](acoustics/geometry.cpp) | 观察点与叶片前后缘坐标变换 |
| | [driver.cpp](acoustics/driver.cpp) | 声学时间状态、边界层表、湍流强度和声学输入读取 |
| `numerics/` | [quadrature.cpp](numerics/quadrature.cpp) | 61 点 Gauss–Kronrod 积分及误差估计 |
| | [backend.cpp](numerics/backend.cpp) | 普通/MKL 后端标识与向量点积 |
| `turbine/aerodynamics/` | [bem.cpp](turbine/aerodynamics/bem.cpp)、[unsteady.cpp](turbine/aerodynamics/unsteady.cpp) | BEM 诱导求解、偏斜修正和非定常翼型状态 |
| `turbine/structure/` | [blade_dynamics.cpp](turbine/structure/blade_dynamics.cpp) | 叶片模态、运动学、广义载荷与加速度 |
| `turbine/coupling/` | [mesh.cpp](turbine/coupling/mesh.cpp) | 气动和结构网格的运动、载荷传递 |
| | [rotor.cpp](turbine/coupling/rotor.cpp) | 叶素状态、气动力和转子装配 |
| | [solver.cpp](turbine/coupling/solver.cpp) | 广义 α 时间推进及结构迭代 |
| `turbine/io/` | [case_input.cpp](turbine/io/case_input.cpp) | 输入解析、配置检查、翼型表及稳态风场 |
| `interfaces/` | [c_api.cpp](interfaces/c_api.cpp) | 对外 C ABI、错误与调用状态 |
| `turbine/simulation/` | [simulation.cpp](turbine/simulation/simulation.cpp) | 完整仿真调度、重置和检查点 |
| | [acoustic_adapter.cpp](turbine/simulation/acoustic_adapter.cpp) | 声学配置装配、气动状态到声学节点的转换 |
| | [results.cpp](turbine/simulation/results.cpp) | 按输出需求分配通道、分块聚合声能 |
| | [batch.cpp](turbine/simulation/batch.cpp) | 独立工况线程调度、输出路径冲突检查和逐项错误 |
| | [file_output.cpp](turbine/simulation/file_output.cpp) | 标准文件输出及失败状态 |
| `apps/` | [section_main.cpp](apps/section_main.cpp) | 单截面频谱示例入口 |
| | [turbine_main.cpp](apps/turbine_main.cpp) | 整机命令行参数解析及运行入口 |
| | [batch_main.cpp](apps/batch_main.cpp) | 批量工况命令行入口 |

`apps/turbine_main.cpp` 调用 `run_case()`；`Simulation` 通过适配器、声学驱动及聚合器处理各时间步，结果交给 `ResultSink`。整机求解器负责协调气动、结构与网格传递；声学模型调用数值积分和后端计算。

声学主路径为 `AcousticDriver::step_blocks()` → `AcousticWorkspace::evaluate_blocks()` → `AcousticAggregator::append()`。每个采样时刻先对各节点调用 `prepare_section()`，再按观察点调用 `emit_section()`。各模型的 `prepare_*()` 计算与观察点无关的谱形，`emit_*()` 施加距离和指向性。源码调用链见 [工作流](../aeroacoustics工作流.md)，性能和数值检查见 [P4–P7 说明](../docs/performance-p4-p7.md)。完整快照接口仍保留。

可靠性辅助接口位于 `include/`：[acoustic_levels.hpp](../include/acoustic_levels.hpp) 区分静音与非法声级；[checked_output.hpp](../include/checked_output.hpp) 检查输出生命周期；[lookup_diagnostics.hpp](../include/lookup_diagnostics.hpp) 提供显式、线程局部的查表诊断会话。掩码、报告和严格模式见 [R1–R3 说明](../docs/reliability.md)。

## 构建文件

- 根目录 `CMakeLists.txt`：构建开关、编译标准、MKL 查找和输出位置。
- [CMakeLists.txt](CMakeLists.txt)：声学静态库及动态库的共用源文件清单。
- [turbine/CMakeLists.txt](turbine/CMakeLists.txt)：整机模块库。
- [apps/CMakeLists.txt](apps/CMakeLists.txt)：截面、整机和批量三个可执行程序。
- 根目录 `tests/CMakeLists.txt`：数值对照探针。

公开头文件路径、构建目标名称及可执行文件位置保持兼容。原有 `section_spectrum()`、`snapshot_spectrum()` 和返回独立快照的 `AcousticDriver::step()` 仍可使用；连续计算可使用复用缓冲区的接口。构建及运行命令见根目录 README。

## 整机工作区与边界层采样

`Solver` 持有 `RotorWorkspace` 和 `LoadWorkspace`。`Rotor::evaluate_into()` 复用气动结果与运动映射数组；`structural_acceleration()` 先重建指定状态的运动，再供载荷映射与质量方程共用。所有暂存区都由调用者持有，不在 `Rotor` 的 `const` 方法中隐藏可变缓存。`transfer_into()` 每次清零累加结果，输入与输出不得使用同一数组。

`PreparedBLTable` 持有已验证表格的私有副本；原 `BLTable` 仍可编辑，其 `interpolate()` 保留校验。整机入口通过 `AcousticDriver::is_sample_time()` 与 `first_node()` 决定哪些节点需要插值。测试及限制见 [P1–P3 说明](../docs/coupling-optimization.md)。

## 模型与状态

`TurbineModel` 深复制 `Case` 并提供只读访问；`Rotor`、UA 翼型引用共享这份模型的所有权。`Solver` 和 `Simulation` 的演化状态私有，复制、移动和完整检查点不会依赖外部容器地址。`AcousticDriver` 的 TI 状态通过 `turbulence_state()` 只读访问。

声源通道信息集中在 [mechanisms.hpp](../include/mechanisms.hpp)，模型选择集中在 `SourceSelection`。TNO 替代压力面、吸力面贡献，保留 BPM 分离贡献。库调用和接口迁移见 [S1–S3 说明](../docs/simulation-api.md)。
