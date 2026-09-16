# P4–P7：声学内存、批量计算与结构求解

实现日期：2026-09-16。对比基线为 `fd193c7`，即 S1–S3 完成后的版本。物理模型和支持范围不变。

## P4：TNO 剖面与频带缓存

`prepare_tno()` 对每个节点的每一侧调用一次 `prepare_profile()`，保存 61 个高度上的边界层量，再供该侧全部频率使用。34 个频率下，剖面计算次数由每侧 34 次降为 1 次。积分矩阵按波数行连续写入，求和仍按原来的高度顺序执行。

频带宽度缓存在 `TnoWorkspace`，频率列表变化时重建。剖面只在当前节点、当前侧的频率循环中复用，每次采样都依据最新边界层重算。直接调用 `spl_integrate()` 仍会准备自身的剖面；没有跨物理状态使用过期缓存。

这项改动不改变 61 点积分规则，也不改变 TNO 外缘速度比约定。整机截面路径仍按参考实现使用单位速度比；直接 `tblte_tno()` 使用传入的速度比。E9 仍是独立待办。

## P5：按需聚合与观察点分块

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

## P6：独立工况并行

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

## P7：求解选项与收敛诊断

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

## 验证与复现

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

## 本机计时

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
