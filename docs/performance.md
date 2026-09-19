# 性能与求解选项

[分块、并行与求解器](#p4-p7) · [结构工作区与边界层](#p1-p3) · [初次优化记录](#initial)

本文集中说明性能接口、缓存边界和历史实测。各阶段保留各自的基线、计时口径与报告，速度结论仅适用于对应工况和机器。工程控制/传播增加了计算量，尚未开展单独的性能优化对比。

<a id="p4-p7"></a>

## 观察点分块、独立工况并行与结构求解（P4–P7）

实现日期：2026-09-16。对比基线为 `fd193c7`，即 S1–S3 完成后的版本。物理模型和支持范围不变。

### P4：TNO 剖面与频带缓存

`prepare_tno()` 对每个节点的每一侧调用一次 `prepare_profile()`，保存 61 个高度上的边界层量，再供该侧全部频率使用。34 个频率下，剖面计算次数由每侧 34 次降为 1 次。积分矩阵按波数行连续写入，求和仍按原来的高度顺序执行。

频带宽度缓存在 `TnoWorkspace`，频率列表变化时重建。剖面只在当前节点、当前侧的频率循环中复用，每次采样都依据最新边界层重算。直接调用 `spl_integrate()` 仍会准备自身的剖面；没有跨物理状态使用过期缓存。

这项改动不改变 61 点积分规则，也不改变 TNO 外缘速度比约定。整机截面路径仍按参考实现使用单位速度比；直接 `tblte_tno()` 使用传入的速度比。E9 仍是独立待办。

### P5：按需聚合与观察点分块

`OutputLayout` 和 `AcousticAggregator` 按 `NrOutFile` 生成标签、分配并累计结果。只请求第 1 类输出时，仍计算所有启用机制对总声能的贡献；第 2–4 类输出数组为空。

整机路径为：

```text
Simulation::next()
  AcousticAggregator::begin()                 # 仅在采样时清零
  AcousticDriver::step_blocks()
    AcousticWorkspace::evaluate_blocks()
      prepare_section()                       # 每节点一次，跨观察点共享
      observe() + emit_section()              # 当前观察点块
      callback -> AcousticAggregator::append()
    TurbulenceState::update()                  # 仍在每个结构步更新
  AcousticAggregator::finish()
  ResultSink::write()
```

默认每块一个观察点，可通过 `RunOptions::observer_block_size` 或命令行 `--observer-block-size=N` 调整。回调按观察点顺序执行，借用的块只在回调期间有效；不要保存引用或从回调重入同一工作区。聚合器检查块的顺序、尺寸和完整性。

`evaluate()`、`step_view()`、`step()`、`snapshot_spectrum()` 继续提供完整快照。直接调用这些接口仍需要容纳全部观察点；分块接口才限制快照峰值。

官方算例每个观察点有 `66 × 7 × 34` 个频谱数值。默认块的频谱有效载荷为 125,664 字节；完整两观察点快照为 251,328 字节。16 个观察点时，完整快照为 2,010,624 字节，分块仍为 125,664 字节。以上只统计 double 数值，未包括容器、声源缓存、积分工作区或输出数组，也不是进程 RSS。

官方两观察点、`NrOutFile=1` 时，聚合数组由原来的 726 个 double 减少为 2 个，即 5,808 字节降为 16 字节。最终请求的输出仍随观察点数量增加。

### P6：独立工况并行

新增 `turbine::run_cases()` 和 `aeroacoustics_batch`。线程领取完整工况，各自创建模型、求解器、TI/UA 状态、查表诊断和 TNO 工作区。返回数组与任务输入顺序一致。一个工况失败时记录错误，其余工况继续完成；命令行只要有一个工况失败便返回非零退出码。

```sh
./build/aeroacoustics_batch 2 case8/IEA_LB_RWT-AeroAcoustics.fst results/wind8 case9/IEA_LB_RWT-AeroAcoustics.fst results/wind9
```

```cpp
#include "turbine/batch.hpp"

turbine::RunOptions options;
options.duration = 20.;
auto outcomes = turbine::run_cases({
    {"case8/IEA_LB_RWT-AeroAcoustics.fst", "results/wind8", options},
    {"case9/IEA_LB_RWT-AeroAcoustics.fst", "results/wind9", options}
}, 2);
```

工作线程数必须为正，实际线程数不超过任务数。启动前检查输出目录重复、嵌套和包含主输入文件的情况；解析已有目录的符号链接或 Windows junction，并处理 Windows 路径大小写。调用方仍须使用专用输出目录，不能在运行中改动输入、输出目录或链接。

MKL 构建在批量工作线程内用 `mkl_set_num_threads_local(1)` 限制内部线程，结束时恢复原线程设置，避免外层工况并行再叠加 MKL 并行。普通版本依赖 C++ 标准线程库，无 OpenMP 依赖。

本轮只并行独立工况。单个工况的时间步、节点和观察点仍串行，保证原求和和状态推进顺序。单工况不会因为选择较大的批量线程数而加速；线程数还应考虑每份模型和工作区所占内存。

### P7：求解选项与收敛诊断

默认 `SolverMode::reference` 保留固定扰动、每轮重建 Jacobian 及修正量判停，与本次改动前的数值路径一致。`SolverOptions` 可设置：

| 字段 | 默认值 | 含义 |
| --- | --- | --- |
| `max_iterations` | 12 | 每个时间步的 Newton 迭代上限 |
| `perturbation` | 1e-4 | 参考模式的加速度差分量；scaled 模式的尺度系数 |
| `correction_tolerance` | 1e-9 | 三片叶片最大修正量范数的阈值 |
| `acceleration_scale` | `{1,1,1}` | 各模态加速度的参考尺度 |
| `residual_absolute` / `residual_relative` | 1e-9 / 1e-9 | 残差的绝对和相对容差 |
| `reuse_jacobian` | true | 仅 scaled 模式使用，限同一时间步内复用 |

`--solver=scaled` 或 `RunOptions::solver.mode = SolverMode::scaled` 启用尺度化模式。第 j 个扰动量为 `perturbation × max(acceleration_scale[j], abs(current_acc[j]))`。停止时要求修正量达标，并同时要求每个分量的以下比值不大于 1：

```text
abs(f[i] - current_acc[i]) /
  (residual_absolute + residual_relative * max(acceleration_scale[i], abs(f[i]), abs(current_acc[i])))
```

Jacobian 只在当前时间步内复用。如果残差没有至少下降 10%，或已经连续复用两次，则重新差分构造；下一时间步一定重建。`reuse_jacobian=false` 可单独核对每轮重建路径。没有引入解析/自动微分 Jacobian，也没有把小型 3×3 系统移交给 MKL。

`Solver::diagnostics()` 返回累计迭代数、加速度计算次数、Jacobian 构造次数、单步最大迭代数，以及最后一轮的残差和修正量。残差是最后一轮线性化点上的值；不是额外重算修正后状态所得。`run.json` 保存模式、累计计数和最后残差。失败消息包含时间、叶片、迭代与范数；整轮达到上限时范数为三片叶片最大值。对状态、加速度、Jacobian、修正量及尺度执行有限值检查。

求解选项、计数和演化状态随复制及检查点保存；`reset()` 保留选项并清零状态与计数。默认值仍不由 `.fst` 的 `ConvTol`、`MaxConvIter` 或 `DT_UJac` 控制。

### 验证与复现

- 与 `fd193c7` 比较：官方 8 m/s、9 m/s、`BLMod=2` 各 20 秒，以及 TNO 0.1 秒，默认模式的数值文件、掩码和查表报告逐字节一致。新增结构诊断字段单独核对。
- 普通和 MKL 后端执行原 Fortran 声学公式、8/9 m/s 整机对照。scaled 模式也执行两个完整整机对照，沿用原容差。
- BPM/TNO 的完整快照与 1/4/13/20 观察点块逐项一致；四类输出需求的共同结果一致。声学库的 TI 方法 1/2 另比较分块与完整路径，覆盖采样、非采样和统计缓冲区填满后的状态。
- 结构测试包含 8/9 m/s 与 0.003125/0.00625/0.0125 s 三种步长、同一步内复用/重建 Jacobian、非零演化状态的检查点重放、重置和不收敛诊断。这不是任意初始挠度或新机型验证。
- 串行与多线程的 BPM/TNO 批量输出逐字节一致；测试独立失败和目录冲突。CI 使用 GCC ThreadSanitizer 检查普通版本的并行路径，不检测 MKL 动态库内部。

```sh
./build/performance_api_probe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/performance-api
python tests/check_performance.py build build/performance-cli
python tests/check_refactor.py /baseline/build/aeroacoustics_turbine build/aeroacoustics_turbine build/before-after --baseline-commit fd193c7f9fb6b467e5f9f06c25e1b187f327691c
python tests/benchmark_performance.py /baseline/build build build/timing --rounds 5
```

测试目录使用新目录。计时应在没有其他构建和测试运行时进行；比较程序使用同一编译器、Release 配置和 MKL 设置。声学内核统计预热后的 C++ 循环耗时，整机与批量统计包含启动和文件输出的外部墙钟时间。

### 本机计时

Windows / GCC 16.2.0，Release，未启用 MKL。新旧版本交替运行五轮，以下为中位数。整机含文件输出；TNO 行为 30 节点、34 频率的声学工作区测试。不同阶段的绝对耗时不能与旧文档直接拼接比较。

| 测试 | 优化前 / 第一种设置 | 优化后 / 第二种设置 | 耗时减少 |
| --- | ---: | ---: | ---: |
| TNO，2 个观察点，3 次频谱计算 | 1.143 s | 1.119 s | 2.1% |
| TNO，16 个观察点，3 次频谱计算 | 1.149 s | 1.132 s | 1.5% |
| 20 s 整机，仅总声压级 | 1.505 s | 1.499 s | 0.4% |
| 20 s 整机，全部输出 | 1.696 s | 1.647 s | 2.8% |
| 20 s 整机，reference → scaled | 1.653 s | 1.528 s | 7.6% |
| 4 个完整工况，1 → 4 个工作线程 | 6.634 s | 2.852 s | 57.0% |

TNO 的变化较小；两组 BPM 频谱测试也有百分之几的波动，不据此宣称普遍提速。默认整机模式的主要收益是减少频谱暂存与无用聚合，耗时接近原版。独立工况并行的吞吐提升更明显。

官方 20 s 工况中，scaled 模式把结构加速度计算从 115,200 次降为 57,984 次，Jacobian 构造从 28,800 次降为 9,696 次；总迭代轮数从 9,600 增为 9,632。整机还包含 BEM、UA、声学和文件输出，因此整体耗时不会按结构调用数的比例下降。

原始五轮计时、可执行文件哈希、公式和整机误差见 [validation-performance.json](validation-performance.json)。结果反映当前机型、工况和硬件，不是跨机器加速保证。

<a id="p1-p3"></a>

## 结构运动、缓冲区与边界层采样（P1–P3）

2026-09-14 完成。优化基线为 [`e5c57ab`](https://github.com/MegrezLi/Aeroacoustics_cpp/commit/e5c57aba793b889e3adf6e8396c9004b5d56ca99)，包含 R1–R3 可靠性修改。本次没有改变气动或声学公式、结构自由度、时间积分系数及 Newton 收敛条件。

### 改动

| 项目 | 修改后的行为 | 主要文件 |
| --- | --- | --- |
| P1 | 叶片基底在气动求值和结构迭代入口分别按叶片准备；一次结构状态求值只建立一组节点运动，载荷映射和加速度装配共用 | `blade_dynamics.cpp`、`rotor.cpp`、`solver.cpp` |
| P2 | `Solver` 持有 `RotorWorkspace`、`LoadWorkspace`；气动结果、结构运动、映射输入和结果数组重复使用 | `rotor.hpp`、`solver.hpp`、`mesh.cpp` |
| P3 | 冻结并验证边界层表；仅在声学采样时，对实际发声节点插值；采样判断由 `AcousticDriver` 统一提供 | `aeroacoustics.hpp`、`driver.cpp`、`turbine_main.cpp` |

#### P1：共享范围明确到一次状态求值

`Rotor::structural_acceleration()` 依次调用 `motions_into()`、`structural_loads_into()` 和接收节点运动的 `acceleration()` 重载。后两步使用同一组运动。每个 Jacobian 扰动都会重建运动，避免把基础状态或上一个扰动的结果带入当前计算。

本算例每片叶片有 17 个结构积分节点，加根、尖共 19 个运动节点。原来一次结构求值在载荷映射中计算 17 次运动，在加速度装配中计算 18 次；现在批量计算 19 次。基底也不再随每个节点重复求三角函数。

气动仍由预测状态计算一次；结构 Newton 校正期间气动分布载荷固定，映射力臂随当前结构状态更新。声学仍读取该步保存的预测气动状态。三片叶片的基础状态和扰动现在逐片处理；现有固定塔架模型中叶片之间没有结构耦合。

#### P2：复用容量，清零累加量

新增 `Rotor::evaluate_into()`、`MotionMap::transfer_into()` 和 `LoadMap::transfer_into()`。原有返回值接口保留，内部调用新实现；需要连续计算时使用工作区接口。

每次映射均重置结果，气动求值先清零平均速度，再写入各节点字段。输入输出使用同一向量时，映射接口报错。工作区由调用者持有，`Rotor` 没有新增隐藏的可变缓存；不同求解器实例使用独立暂存区。数组长度变化时允许重新分配，不保证所有未来网格和异常路径都没有分配。

风场节点查询的错误字符串改为仅在异常分支构造，避免部分编译器在每次有效查询时也分配临时字符串；地面以下的幂律风查询仍然报错。

`turbine_workspace_probe` 对本算例预热后连续 50 个完整 `Solver::step()` 计数，观察到 **0 次 `operator new/new[]` 调用**。该计数不包括初始化、文件输出或异常诊断，也不代表程序完全不占用堆内存。

#### P3：冻结表格，按声学需求查询

`PreparedBLTable` 保存表格的私有副本，构造时检查坐标轴、数组形状及有限值；后续查询只检查查询参数并执行范围诊断和插值。原 `BLTable` 的公开数组仍可编辑，旧 `interpolate()` 保留每次检查。修改原表不会影响已准备的副本，需要更新时重新构造。

`is_sample_time()` 和 `first_node()` 给整机入口提供与声学驱动一致的采样条件。默认 20 秒、`DT=0.00625`、`DT_AA=0.1`、`BldPrcnt=70` 下，如果启用表格边界层，调用次数由 `3201 × 90 = 288090` 降为 `201 × 66 = 13266`，减少约 **95.4%**。这是按循环与配置计算的插值次数，不能直接解释为整机加速比例；默认经验边界层工况没有这部分插值。

所有节点的速度、入流和湍流强度仍按每个时间步更新，声学采样使用更新前的 TI。边界层超界报告现在只统计实际执行的声学查询；严格模式也只在这些查询发生时检查范围。整机初始化仍读取并验证所需表文件。BEM 和 UA 的查表诊断范围不变。

### 实测

Windows、GCC 16.2.0、Release `-O3 -DNDEBUG`，未启用 MKL。旧程序和新程序分别预热一次，再交替运行五轮，取内部墙钟时间的中位数，包含输入读取和文件输出。三个计时工况均运行 20 秒、3200 步，输出 201 个声学时刻。

| 工况 | 优化前 | 优化后 | 耗时减少 |
| --- | ---: | ---: | ---: |
| 官方 8 m/s | 1.985 s | 1.627 s | 18.1% |
| 9 m/s 扰动 | 2.005 s | 1.570 s | 21.7% |
| 8 m/s，`BLMod=2` | 2.188 s | 1.749 s | 20.0% |

这组数据衡量 P1–P3 的整体收益，没有单独归因到每一项，也没有测量本次 MKL 加速比。结果受机器负载、编译器及输出规模影响。逐轮时间、程序和输出 SHA256 见 [验证数据](validation-coupling.json)。之前一轮优化的数据仍见 [原优化记录](performance.md#initial)，两轮基线不同。

### 数值与接口检查

- 上述三个工况，以及 `TBLTEMod=2` 的 0.1 秒短工况：每个工况的 `dynamics.csv`、四份声学 `.out` 和四份 `.mask` 与优化前逐字节一致。短 TNO 工况只用于回归，没有作为计时基准。
- 独立结构探针的 5700 行、网格探针的 4700 行，以及 2 秒求解探针的全部节点气动力和状态，与旧程序逐字节一致。
- 普通和 MKL 版本均通过官方 8 m/s、9 m/s 整机 Fortran 对照；每个工况核对 101 个输入文件和 145926 个声学值，原容差保持不变。
- 两个后端均通过 204 个声学配置、56576 个数值的 Fortran 公式对照，包含 16 个 TNO 配置，最大绝对误差约 `4.55×10⁻¹³ dB`。
- 新增探针检查重复映射清零、刚体平移、恒定线载荷守恒、输入输出别名拒绝、只读表与原表分离、双线性插值及端点钳制、严格查表策略、采样时间，以及两种 TI 方法的逐步更新和采样前状态。
- R1–R3 的单元及命令行异常测试通过。GitHub Actions 增加了整机工作区探针。

这些结果验证实现一致性，不替代风机噪声的实测验证。

### 复现

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build -j 4
./build/turbine_workspace_probe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst
python tests/benchmark_coupling.py /path/to/baseline/aeroacoustics_turbine build/aeroacoustics_turbine build/coupling-benchmark --baseline-commit e5c57aba793b889e3adf6e8396c9004b5d56ca99
```

基准程序应从所指定提交单独构建，并与候选程序使用相同编译器和优化选项。输出目录须不存在，脚本会在该目录复制工况并修改风速或模型开关。Windows 给程序路径加 `.exe`；运行前配置编译器运行库路径。Python 只负责任务调度和结果比较。

<a id="initial"></a>

## 声源准备与发射分离：初次优化记录

本次优化针对结构载荷重复映射、声源重复计算、临时数组分配、固定配置重复解析，以及经验模型中沿用的 Fortran 式组织方式。基线为 `ad44bfd`，物理模型、输入文件、时间网格和收敛阈值保持不变。

### 修改内容

| 问题 | 修改 | 主要代码 |
| --- | --- | --- |
| 数值 Jacobian 扰动一片叶片时，重新映射三片叶片的载荷 | 新增单叶片映射函数；每轮 Newton 迭代由 30 次单叶片映射减少为 12 次 | [rotor.cpp](../src/turbine/coupling/rotor.cpp)、[solver.cpp](../src/turbine/coupling/solver.cpp) |
| 每个观察点重复计算同一节点的边界层、声源谱形和 TNO 积分 | 分为 `prepare_*()` 与 `emit_*()`；先按节点准备声源，再按观察点施加距离和指向性 | [source_models.hpp](../src/acoustics/source_models.hpp)、[spectrum.cpp](../src/acoustics/spectrum.cpp)、[workspace.cpp](../src/acoustics/workspace.cpp) |
| 时间循环、频谱装配和 TNO 积分反复创建或复制数组 | 工作区保存节点、频谱和积分数组；主程序复用四类输出数组，并通过只读快照指针读取结果 | [driver.cpp](../src/acoustics/driver.cpp)、[tno.cpp](../src/acoustics/tno.cpp)、[turbine_main.cpp](../src/apps/turbine_main.cpp) |
| 节点和输出循环反复解析固定输入，重复校验参数和计算 A 计权 | 初始化时缓存 `SkewOptions`、`NrOutFile`、已校验声学参数和 A 计权；动态截面输入仍逐次检查 | [rotor.cpp](../src/turbine/coupling/rotor.cpp)、[workspace.cpp](../src/acoustics/workspace.cpp) |
| 经验公式集中在一个大文件，参数列表长，变量声明和数组下标沿用 Fortran 风格 | 按边界层、BPM 尾缘、其他 BPM 声源、入流、谱形和指向性拆分；使用具名数据类型、局部变量及从零开始的循环 | [源码索引](development.md#source)、[kernel_compat.cpp](../src/acoustics/kernel_compat.cpp) |

BPM 的共享边界层在每个节点、每个声学采样时刻只计算一次。TNO 两侧积分也只准备一次，观察点变化只影响后续距离和指向性计算。缓存不会跨时间保留过期声源：新的速度、攻角、湍流强度和边界层输入会重新生成当前谱形。

经验公式的系数、分段条件、钝度频谱逐频带累积归一化和声级下限沿用原实现。此次没有增加噪声模型或改变气动与结构的推进次序。

### 性能结果

环境：Windows 11、Intel Core i7-11800H、GCC 16.2.0、Release（`-O3 -DNDEBUG`），普通 C++ 后端，未启用 MKL。基线和优化版先预热，再交替运行五轮，以下取中位数。原始样本、可执行文件 SHA256 和输出差异保存在 [validation-optimization.json](validation-optimization.json)。

#### 完整 20 秒算例

| 工况 | 优化前 | 优化后 | 加速比 | 耗时减少 |
| --- | ---: | ---: | ---: | ---: |
| 官方 8 m/s | 3.536 s | 1.824 s | 1.94× | 48.4% |
| 9 m/s 扰动 | 3.497 s | 1.794 s | 1.95× | 48.7% |

计时来自 `run.json`，覆盖程序内的初始化、3200 步推进、201 次声学采样及输出处理，不包含外部启动进程的开销。没有缩短物理时长或减少输出。

#### 声学局部测试

每次快照使用 30 个节点和 34 个频率；BPM 测量批次包含 20 次快照，TNO 包含 2 次，表内折算为一次快照。基线调用原 `snapshot_spectrum()`，优化版复用 `AcousticWorkspace`，因此加速同时包含声源共享和缓冲区复用的收益。

| 模型 | 观察点数 | 优化前 | 优化后 | 加速比 |
| --- | ---: | ---: | ---: | ---: |
| BPM | 2 | 2.013 ms | 0.659 ms | 3.06× |
| BPM | 16 | 16.508 ms | 1.169 ms | 14.13× |
| TNO | 2 | 739.289 ms | 361.894 ms | 2.04× |
| TNO | 16 | 6078.342 ms | 364.610 ms | 16.67× |

多观察点共享计算时收益更明显。以上是特定测试输入的局部结果，不能直接作为整机加速比。官方整机默认使用 BPM，不启用 TNO；此次数据也不是 MKL 与普通后端的速度比较。

### 数值检查

- **Fortran 公式对照**：普通和 MKL 后端分别通过 56,576 个值、204 组配置的对照，包含 16 组 TNO 配置。最大绝对差异约 `4.55×10⁻¹³ dB`，阈值为 `2×10⁻⁸ dB`。见 [普通后端报告](validation-fortran-portable.json) 和 [MKL 报告](validation-fortran-mkl.json)。
- **整机 Fortran 对照**：官方 8 m/s 和 9 m/s 扰动各检查 145,926 个声学值及 101 个输入文件，均满足原有阈值。见 [官方工况报告](validation-full-case.json) 和 [扰动工况报告](validation-wind9.json)。
- **优化前后整机对照**：两个工况的 `dynamics.csv` 均逐字节一致；9 m/s 的四份声学文件也逐字节一致。8 m/s 仅频带输出文件出现约 `1×10⁻¹⁰ dB` 的末位变化，其余三份一致。
- **工作区复用**：36 组配置与状态组合改变节点数、观察点数、攻角、TI、边界层和 A 计权，检查复用输出与独立截面调用，并与旧 C++ 库比较 21,420 个值；有限值最大差异约 `5.68×10⁻¹⁴ dB`。普通和 MKL 后端均通过复用检查。
- **声学驱动状态**：两个 TI 计算分支、边界层表、展向选点、厚度提取和积分辅助功能与旧 C++ 输出对照通过；驱动频谱最大差异约 `1.14×10⁻¹³ dB`。

公式对照误差与整机误差含义不同。整机仍包含 C++ 和原 Fortran 求根及耦合迭代的数值差异，原有支持范围见 [模块说明](development.md#scope)。

### 接口与复现

`AcousticWorkspace::evaluate()` 返回内部快照的只读引用，`AcousticDriver::step_view()` 在采样时刻返回只读指针，其他时刻返回空指针。调用者须在下一次相应调用前使用完数据；需要持久保存时自行复制。各求解实例分别持有工作区，同一工作区不应被多个线程同时写入。

原 `AcousticDriver::step()` 继续返回独立快照；`section_spectrum()`、`snapshot_spectrum()`、标量模型函数和 C ABI 保留。应用使用新 C++ 头文件后应重新编译。整机仍可仅用 C++ 编译和运行；Python 脚本仅负责测试调度和数值文件比较。

构建与核心验证：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build -j 4
./build/acoustic_workspace_probe verify build/workspace.csv
python tests/compare_fortran.py build/libaeroacoustics_shared.so --fortran-compiler gfortran --report build/kernels.json
python tests/run_standalone.py build/aeroacoustics_turbine build/official-check --report build/full-case.json
python tests/run_perturbation.py build/aeroacoustics_turbine build/wind9 --report build/wind9.json
```

在同一机器上构建优化前基线，并将新探针以兼容模式链接至旧库。下面使用单配置生成器、GCC/Clang：

```sh
git worktree add --detach build-baseline-src ad44bfd
cmake -S build-baseline-src -B build-baseline -DCMAKE_BUILD_TYPE=Release
cmake --build build-baseline -j 4
c++ -std=c++17 -O3 -DNDEBUG -DAERO_BENCHMARK_LEGACY -Ibuild-baseline-src/include tests/acoustic_workspace_probe.cpp build-baseline/libaeroacoustics.a -o build-baseline/acoustic_workspace_probe
python tests/benchmark_optimization.py build-baseline build --wind9-case build/wind9/inputs/IEA_LB_RWT-AeroAcoustics.fst --repeats 5 --baseline-commit ad44bfd --report build/optimization.json
```

Windows / MinGW 先运行 `. ./tools/use_gcc.ps1`，使用 `-G 'MinGW Makefiles'`，将 `.so` 改为 `libaeroacoustics_shared.dll` 并为可执行文件加 `.exe`；手动编译探针时将 `c++` 换为 `g++`。9 m/s 输入由前面的 `run_perturbation.py` 生成。复现实测值时应保持编译选项一致，并避免同时执行其他高负载任务。
