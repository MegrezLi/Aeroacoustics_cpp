# S1–S3：仿真接口、声源配置与状态管理

本轮基线为 [`2d7fda5`](https://github.com/MegrezLi/Aeroacoustics_cpp/commit/2d7fda532b389371a854c49daad9e96282c938b1)。计算公式、时间步长、Newton 收敛条件和七类标准声学通道保持原有定义。

## S1：仿真与输出分开

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

## S2：机制与模型分别描述

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

## S3：模型所有权与完整状态

`TurbineModel` 深复制传入的 `Case`，向外只提供 const 访问。`Rotor` 共享这份模型；其非定常翼型对象通过共享所有权引用其中的翼型，不依赖调用者容器的地址。单独构造 `UnsteadyAirfoil(const Airfoil&, ...)` 时，该对象复制并持有翼型数据。

`Solver` 的结构状态、气动结果、时间、加速度和转子改为私有成员。常用读取方式改为：

| 原调用 | 当前调用 |
| --- | --- |
| `solver.state` | `solver.state()` |
| `solver.aerodynamic` | `solver.aerodynamic()` |
| `solver.time`、`solver.step_number` | `solver.time()`、`solver.step_number()` |
| `solver.rotor.structure` | `solver.rotor().structure()` |
| `driver.state.values` | `driver.turbulence_state().values` |

这些访问器返回值或只读引用，不能单独改写时间或位移，绕过相应的气动和历史状态更新。需要不同初始模型时构造新的 `Case`/`TurbineModel`。本轮没有增加任意初始挠度等物理能力。

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

`run.json` 的物理时间、步数和采样数保持不变；`elapsed_seconds` 现在统计从 `Simulation` 初始化到最后一帧处理完毕的墙钟时间，不含此前的 `Case` 文件解析及最后的输出关闭。跨版本性能比较应使用外部墙钟计时。本轮以结构和接口验证为主，不宣称额外加速。

## 验证

- 官方 8 m/s、9 m/s 和 `BLMod=2` 各运行 20 秒；另运行 0.1 秒 TNO 工况。结构、声学、掩码和查表报告与基线逐字节一致。除墙钟时间外，运行元数据一致。
- 普通与 MKL 版本均执行 Fortran 公式和两个整机工况对照，保持原容差。
- 检查 216 组声源开关组合、标准通道顺序及 TNO 替代关系。
- 库接口直接运行与分段运行结果完全一致；模型销毁、外部数据修改、复制、移动、检查点恢复、重置和诊断会话隔离均有测试。
- 将测试输入副本移走后，已构造的仿真仍能重置并重现结果；官方输入不作修改。
- P1–P3 的工作区检查继续通过，预热后的 50 个 `Solver::step()` 无 `operator new/new[]` 调用；R1–R3 异常检查保留。

报告见 [validation-structure.json](validation-structure.json)。主要复现命令：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build -j 4
./build/mechanism_probe
./build/simulation_api_probe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/api-check
python tests/check_refactor.py /path/to/baseline/aeroacoustics_turbine build/aeroacoustics_turbine build/refactor-check
```

比较基线为本文件开头的提交；两个比较程序应使用相同编译器和构建选项。测试输出目录使用新目录，Windows 给程序名加 `.exe`。Python 只负责测试调度和文件比较。

## P4–P7 接口补充

`RunOptions` 新增 `observer_block_size`（默认 1）和 `solver`（`SolverOptions`）。`Simulation` 的声学路径改用观察点分块回调，`AcousticResult::power` 中未由 `NrOutFile` 请求的类别为空。需要完整四类结果时设置 `NrOutFile=4`；已有完整快照 API 不变。

`Solver::diagnostics()` 与 `RunSummary::structure` 提供求解统计。默认参考模式保持原数值路径，尺度化模式须明确选择。`run.json` 新增 `structural_*` 字段。`SolverOptions` 随检查点恢复，`reset()` 保留选项。

`turbine/batch.hpp` 提供独立工况的 `run_cases()`。线程、路径限制、借用块生命周期、完整参数与验证见 [P4–P7 说明](performance-p4-p7.md)。
