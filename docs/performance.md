# 性能与求解选项

<a id="p4-p7"></a>

## 观察点分块、工况并行与结构求解

### TNO 剖面与频带缓存

`prepare_tno()` 对每个节点的每一侧调用一次 `prepare_profile()`，保存 61 个高度上的边界层量，再供该侧全部频率使用。34 个频率下，剖面计算次数由每侧 34 次降为 1 次。积分矩阵按波数行连续写入，求和仍按原来的高度顺序执行。

频带宽度缓存在 `TnoWorkspace`，频率列表变化时重建。剖面只在当前节点、当前侧的频率循环中复用，每次采样都依据最新边界层重算。直接调用 `spl_integrate()` 仍会准备自身的剖面；没有跨物理状态使用过期缓存。

TNO 使用 61 点积分。整机默认 `reference` 模式采用单位外缘速度比；`input` 模式使用边界层输入速度比，配置见 [尾缘模型](engineering.md#trailing-edge)。

### 按需聚合与观察点分块

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

### 独立工况并行

`turbine::run_cases()` 和 `aeroacoustics_batch`。线程领取完整工况，各自创建模型、求解器、TI/UA 状态、查表诊断和 TNO 工作区。返回数组与任务输入顺序一致。一个工况失败时记录错误，其余工况继续完成；命令行只要有一个工况失败便返回非零退出码。

先准备 `case8/` 和 `case9/` 两份独立输入，并创建父目录 `results/`：

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

批量接口并行独立工况。单个工况的时间步、节点和观察点仍串行，保证原求和和状态推进顺序。单工况不会因为选择较大的批量线程数而加速；线程数还应考虑每份模型和工作区所占内存。

### 求解选项与收敛诊断

默认 `SolverMode::reference` 使用固定扰动、每轮重建 Jacobian 及修正量判停。`SolverOptions` 可设置：

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

Jacobian 只在当前时间步内复用。如果残差没有至少下降 10%，或已经连续复用两次，则重新差分构造；下一时间步一定重建。`reuse_jacobian=false` 可单独核对每轮重建路径。

`Solver::diagnostics()` 返回累计迭代数、加速度计算次数、Jacobian 构造次数、单步最大迭代数，以及最后一轮的残差和修正量。残差是最后一轮线性化点上的值；不是额外重算修正后状态所得。`run.json` 保存模式、累计计数和最后残差。失败消息包含时间、叶片、迭代与范数；整轮达到上限时范数为三片叶片最大值。对状态、加速度、Jacobian、修正量及尺度执行有限值检查。

求解选项、计数和演化状态随复制及检查点保存；`reset()` 保留选项并清零状态与计数。默认值仍不由 `.fst` 的 `ConvTol`、`MaxConvIter` 或 `DT_UJac` 控制。

<a id="p1-p3"></a>

## 缓存与内存

| 对象 | 复用范围与边界 |
| --- | --- |
| 结构运动、气动网格 | 一次状态求值共享；状态改变后重算，Newton 试探状态不能复用上一轮结果 |
| `PreparedBLTable` | 校验一次并持有独立表格；插值仍使用当前攻角、Re 和弦长 |
| `SurfaceSet` | 元数据只读一次，按翼型编号索引查找；表面数据顺序保持输入顺序 |
| `TnoWorkspace` | 首次有效 TNO 求值才分配积分数组；BPM/Howe 不承担这部分内存 |
| 观察点分块 | 完整块和尾块分别复用；固定配置预热后不再反复分配内部频谱 |
| `ArrivalHistory` | 连续时间/声能数组，插值直接累加到目标缓冲区；查询游标属于调用者 |
| 室外传播 | 同一节点分别准备前缘、尾缘两条路径，供对应声源共享；不跨时间缓存 |
| TI 方法 2 | 滚动均值和中心平方和，每个窗口长度完整重算一次，摊销更新为 O(1) |

分块大小为 B 且存在尾块时，最多保留 B 个完整块位置和 B−1 个尾块位置。连续历史仍随时长增长，扩容会复制数据；`MaxValues` 继续限制容量。累计分配字节不等于峰值内存。

工作区归各个仿真对象所有，复制和检查点保存独立状态。TI 滚动中心矩使用 `long double`，与逐次全窗口求和允许舍入差异。C++ 对象布局变化后，调用方需重新编译。

## 验证与性能比较

```sh
ctest --test-dir build -L unit --output-on-failure
python tests/benchmark_optimizations.py /baseline/optimization_benchmark /candidate/optimization_benchmark build/paired-benchmark --rounds 5
python tests/check_refactor.py /baseline/aeroacoustics_turbine /candidate/aeroacoustics_turbine build/refactor-check --baseline-commit BASELINE_SHA
```

基准程序分别用对应版本的头文件和静态库编译，两边使用相同编译器、Release 选项与 MKL 设置。脚本预热后交替运行，记录中位数、分配次数、累计字节、数值校验和与程序哈希。`check_refactor.py` 比较默认 8 m/s、9 m/s、表格边界层和短时 TNO 工况的数值文件。

`check_refactor.py` 加 `--rounds 5` 可预热并交替计时；工作区、输出数量、求解器和批量运行的专项计时使用 `tests/benchmark_performance.py BEFORE_BUILD AFTER_BUILD OUTPUT --baseline-commit BASELINE_SHA`。

比较时固定硬件、编译选项和工况，使用脚本输出的配对计时。历史原始数据：[缓存与分配](validation-optimization-five.json)、[分块与求解器](validation-performance.json)。
