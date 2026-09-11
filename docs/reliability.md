# R1–R3：声级、输出和查表诊断

本次处理非有限声级、零声能输出语义、文件写入检查和查表范围诊断。物理公式及默认端点插值策略保持不变。

## R1：区分静音与计算异常

- `-Inf` 继续表示零声能，参与总量求和时贡献为零。
- `NaN`、`+Inf`、声级转声能时溢出，以及累加得到的非有限/负声能均报错。
- 声源模型产生非法声级时，整机错误信息包含时间、观察点、叶片、原始节点编号、机制、频率和尾缘/入流模型开关。
- 单截面 C++ 接口报告机制与频率；工作区额外报告观察点和节点。独立库用户未指定诊断编号时，编号为 0，另提供工作区内的节点序号。

关键实现：[acoustic_levels.hpp](../include/acoustic_levels.hpp)、[spectrum.cpp](../src/acoustics/spectrum.cpp)、[workspace.cpp](../src/acoustics/workspace.cpp)、[driver.cpp](../src/acoustics/driver.cpp)。异常会使整机返回非零退出码，不再把它当作静音继续输出。

## R2：保持旧文件，同时给出零声能掩码

原有四类 `.out` 保留 Fortran 对照格式，其中零声能仍使用 `0 dB` 占位。每份实际生成的 `.out` 都增加同名 `.mask`，行、时间和通道顺序完全对应：

| mask | 含义 | 如何解释 `.out` |
| --- | --- | --- |
| 0 | 当前聚合结果为零声能 | `0 dB` 是占位值；可来自关闭的机制、未选中的节点或数值上为零的贡献 |
| 1 | 当前聚合结果为正声能 | 声级是有效结果；即使它恰好为 `0 dB`，也不是静音占位 |

掩码描述声能是否为零，不表示预测可信度，也不细分零声能的物理原因。非法值会报错，不会写成 mask=0。每份文件仍有三行表头；`run.json` 记录该编码规则。

[CheckedOutput](../include/checked_output.hpp) 对声学输出、掩码、动力学、查表报告和运行记录统一检查打开、写入、刷新与关闭。错误包含失败文件路径；完成提示只在所有输出关闭成功后出现。

开始生成新输出时先清空原 `run.json`，在数据文件均完成后才写入新运行记录。失败可能留下部分数据文件；调用者必须检查退出码，不能仅以 `.out` 存在判断成功。输入/参数解析阶段失败属于尚未开始输出的新运行，应同样以退出码为准。

## R3：默认汇总诊断，可选严格报错

默认策略 `clamp`：超出表格范围后仍使用端点值，并在成功运行后生成 `lookup_diagnostics.csv`；有越界时终端给出一次汇总警告。

```sh
./build/aeroacoustics_turbine examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/results --lookup-policy=clamp
```

严格策略 `error`：首次超界即失败，错误信息包含文件、坐标轴、输入值、表格范围、时间、叶片、节点和查询阶段。

```sh
./build/aeroacoustics_turbine examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/strict 20 --lookup-policy=error
```

时长参数与策略参数可以交换顺序；Windows 可执行文件加 `.exe`。不新增未经验证的外推公式。严格模式也检查 BEM 试探角度和非定常翼型内部查询，因此严格模式失败并不自动说明最终工作攻角越界。

报告按表格、坐标轴、叶片、节点和阶段分组，包含：

- `calls`：实际超界查询次数；同一步重复查询也计数，一个 BL 查询的两个轴越界分别计数。
- `first_time_s`、`last_time_s`：该组超界查询的最早与最晚输入时刻，不是持续越界时长。
- `min_value`、`max_value`、`lower`、`upper`：超界输入极值和表格有效范围。
- `axis`：`alpha_deg`（度）或 `Re`（无量纲）；翼型攻角先沿用原来的周期归一化，再检查表格范围。
- `stage`：`BEM_trial`、`UA_evaluate`、`UA_advance` 或 `BL_interpolate`。它们用于区分求根试探和不同调用阶段，不能直接当作最终物理状态的超界比例。

精确落在端点属于有效输入，不计越界；没有越界时仍生成仅含表头的报告。失败时使用终端定位信息，部分或空的 CSV 不能作为完整统计。没有改变原来的边界层插值调用频率或 TI 更新时间顺序。

库级入口保持原签名。需要报告时，可使用 [lookup_diagnostics.hpp](../include/lookup_diagnostics.hpp) 的 `LookupReport`、`LookupSession` 和 `LookupLocation`。会话是显式启用、线程局部的，嵌套结束后恢复调用者上下文；未启用会话的库调用保持端点策略。多线程用户应为每个线程建立独立报告，再合并结果。

## 验证与复现

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build -j 4
./build/reliability_probe examples/IEA_LB_RWT-AeroAcoustics/Airfoils/RotorSE_FAST_IEA_landBased_RWT_AeroDyn_Polar_10.dat build/reliability-unit
python tests/check_reliability.py build/aeroacoustics_turbine build/reliability-cli
python tests/run_standalone.py build/aeroacoustics_turbine build/official-check --report build/full-case.json
python tests/run_perturbation.py build/aeroacoustics_turbine build/wind9 --report build/wind9.json
```

CLI 测试要求新的输出目录，以保留每次故障测试的证据。测试中修改的输入仅在该目录的副本内，不修改官方输入。

针对性检查覆盖有效 0 dB 与静音、NaN/+Inf、声能溢出、模型异常上下文、真实查表的端点/超界、严格模式、嵌套诊断会话、输出路径冲突及刷新失败。Linux 还使用 `/dev/full` 检查真实写入失败。普通/MKL 声学检查及整机 Fortran 对照继续保留；此次加强可靠性，不作为新的提速声明。

本机验证：普通/MKL 两个后端的可靠性探针与声学工作区检查通过；两个后端各完成 56,576 个声学值的 Fortran 对照。8 m/s 和 9 m/s 整机对照通过，四份 `.out` 及 `dynamics.csv` 与 `3c24ded` 版本逐字节一致；两个工况都未触发查表超界。CLI 故障测试覆盖五类输出路径冲突、模型非有限声级及查表严格模式。结果见 [validation-reliability.json](validation-reliability.json)。掩码和诊断增加了输出及检查工作量，运行时间应在相同输出配置下比较。
