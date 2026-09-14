# P1–P3：结构运动、缓冲区与边界层采样

2026-09-14 完成。优化基线为 [`e5c57ab`](https://github.com/MegrezLi/Aeroacoustics_cpp/commit/e5c57aba793b889e3adf6e8396c9004b5d56ca99)，包含 R1–R3 可靠性修改。本次没有改变气动或声学公式、结构自由度、时间积分系数及 Newton 收敛条件。

## 改动

| 项目 | 修改后的行为 | 主要文件 |
| --- | --- | --- |
| P1 | 叶片基底在气动求值和结构迭代入口分别按叶片准备；一次结构状态求值只建立一组节点运动，载荷映射和加速度装配共用 | `blade_dynamics.cpp`、`rotor.cpp`、`solver.cpp` |
| P2 | `Solver` 持有 `RotorWorkspace`、`LoadWorkspace`；气动结果、结构运动、映射输入和结果数组重复使用 | `rotor.hpp`、`solver.hpp`、`mesh.cpp` |
| P3 | 冻结并验证边界层表；仅在声学采样时，对实际发声节点插值；采样判断由 `AcousticDriver` 统一提供 | `aeroacoustics.hpp`、`driver.cpp`、`turbine_main.cpp` |

### P1：共享范围明确到一次状态求值

`Rotor::structural_acceleration()` 依次调用 `motions_into()`、`structural_loads_into()` 和接收节点运动的 `acceleration()` 重载。后两步使用同一组运动。每个 Jacobian 扰动都会重建运动，避免把基础状态或上一个扰动的结果带入当前计算。

本算例每片叶片有 17 个结构积分节点，加根、尖共 19 个运动节点。原来一次结构求值在载荷映射中计算 17 次运动，在加速度装配中计算 18 次；现在批量计算 19 次。基底也不再随每个节点重复求三角函数。

气动仍由预测状态计算一次；结构 Newton 校正期间气动分布载荷固定，映射力臂随当前结构状态更新。声学仍读取该步保存的预测气动状态。三片叶片的基础状态和扰动现在逐片处理；现有固定塔架模型中叶片之间没有结构耦合。

### P2：复用容量，清零累加量

新增 `Rotor::evaluate_into()`、`MotionMap::transfer_into()` 和 `LoadMap::transfer_into()`。原有返回值接口保留，内部调用新实现；需要连续计算时使用工作区接口。

每次映射均重置结果，气动求值先清零平均速度，再写入各节点字段。输入输出使用同一向量时，映射接口报错。工作区由调用者持有，`Rotor` 没有新增隐藏的可变缓存；不同求解器实例使用独立暂存区。数组长度变化时允许重新分配，不保证所有未来网格和异常路径都没有分配。

`turbine_workspace_probe` 对本算例预热后连续 50 个完整 `Solver::step()` 计数，观察到 **0 次 `operator new/new[]` 调用**。该计数不包括初始化、文件输出或异常诊断，也不代表程序完全不占用堆内存。

### P3：冻结表格，按声学需求查询

`PreparedBLTable` 保存表格的私有副本，构造时检查坐标轴、数组形状及有限值；后续查询只检查查询参数并执行范围诊断和插值。原 `BLTable` 的公开数组仍可编辑，旧 `interpolate()` 保留每次检查。修改原表不会影响已准备的副本，需要更新时重新构造。

`is_sample_time()` 和 `first_node()` 给整机入口提供与声学驱动一致的采样条件。默认 20 秒、`DT=0.00625`、`DT_AA=0.1`、`BldPrcnt=70` 下，如果启用表格边界层，调用次数由 `3201 × 90 = 288090` 降为 `201 × 66 = 13266`，减少约 **95.4%**。这是按循环与配置计算的插值次数，不能直接解释为整机加速比例；默认经验边界层工况没有这部分插值。

所有节点的速度、入流和湍流强度仍按每个时间步更新，声学采样使用更新前的 TI。边界层超界报告现在只统计实际执行的声学查询；严格模式也只在这些查询发生时检查范围。整机初始化仍读取并验证所需表文件。BEM 和 UA 的查表诊断范围不变。

## 实测

Windows、GCC 16.2.0、Release `-O3 -DNDEBUG`，未启用 MKL。旧程序和新程序分别预热一次，再交替运行五轮，取内部墙钟时间的中位数，包含输入读取和文件输出。三个计时工况均运行 20 秒、3200 步，输出 201 个声学时刻。

| 工况 | 优化前 | 优化后 | 耗时减少 |
| --- | ---: | ---: | ---: |
| 官方 8 m/s | 1.985 s | 1.627 s | 18.1% |
| 9 m/s 扰动 | 2.005 s | 1.570 s | 21.7% |
| 8 m/s，`BLMod=2` | 2.188 s | 1.749 s | 20.0% |

这组数据衡量 P1–P3 的整体收益，没有单独归因到每一项，也没有测量本次 MKL 加速比。结果受机器负载、编译器及输出规模影响。逐轮时间、程序和输出 SHA256 见 [验证数据](validation-coupling.json)。之前一轮优化的数据仍见 [原优化记录](optimization.md)，两轮基线不同。

## 数值与接口检查

- 上述三个工况，以及 `TBLTEMod=2` 的 0.1 秒短工况：每个工况的 `dynamics.csv`、四份声学 `.out` 和四份 `.mask` 与优化前逐字节一致。短 TNO 工况只用于回归，没有作为计时基准。
- 独立结构探针的 5700 行、网格探针的 4700 行，以及 2 秒求解探针的全部节点气动力和状态，与旧程序逐字节一致。
- 普通和 MKL 版本均通过官方 8 m/s、9 m/s 整机 Fortran 对照；每个工况核对 101 个输入文件和 145926 个声学值，原容差保持不变。
- 两个后端均通过 204 个声学配置、56576 个数值的 Fortran 公式对照，包含 16 个 TNO 配置，最大绝对误差约 `4.55×10⁻¹³ dB`。
- 新增探针检查重复映射清零、刚体平移、恒定线载荷守恒、输入输出别名拒绝、只读表与原表分离、双线性插值及端点钳制、严格查表策略、采样时间，以及两种 TI 方法的逐步更新和采样前状态。
- R1–R3 的单元及命令行异常测试通过。GitHub Actions 增加了整机工作区探针。

这些结果验证实现一致性，不替代风机噪声的实测验证。

## 复现

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build -j 4
./build/turbine_workspace_probe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst
python tests/benchmark_coupling.py /path/to/baseline/aeroacoustics_turbine build/aeroacoustics_turbine build/coupling-benchmark --baseline-commit e5c57aba793b889e3adf6e8396c9004b5d56ca99
```

基准程序应从所指定提交单独构建，并与候选程序使用相同编译器和优化选项。输出目录须不存在，脚本会在该目录复制工况并修改风速或模型开关。Windows 给程序路径加 `.exe`；运行前配置编译器运行库路径。Python 只负责任务调度和结果比较。
