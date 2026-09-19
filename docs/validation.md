# 验证与复现

[Fortran 对照](#fortran) · [声级、输出与查表](#reliability) · [接口与并行](#interfaces) · [工程模型](#engineering)

Fortran 对照检验移植一致性；解析测试与步长收敛检验新增模型的实现。两者都不能替代现场声学测量和机型参数标定。所有命令在仓库根目录执行，测试输出目录应使用新目录。

<a id="fortran"></a>

## Fortran 与独立 C++ 对照

所有数值对照均以原 Fortran 为参考。默认 C++ 求解器不使用 `reference/`；测试脚本在运行结束后读取参考结果并比较。

### 随仓库运行回归

先按 README 编译。以下命令在项目根目录执行，需要 Python 3.11+ 与 NumPy。Windows 对可执行文件加 `.exe`。

```sh
python tests/run_standalone.py build/aeroacoustics_turbine build/official-check --report build/full-case.json
python tests/run_perturbation.py build/aeroacoustics_turbine build/wind9-check --report build/wind9.json
```

第一项核对官方 101 个输入文件的 SHA256，再运行完整 20 秒。第二项复制相同输入，只将 InflowWind 的 `HWindSpeed` 改成 9 m/s，并与单独运行的原 Fortran 对比。两项均逐值比较全部四个声学输出文件，共 201 个时刻、145,926 个声学值。

| 输出 | 通道数 | 绝对容差 |
| --- | ---: | ---: |
| 总声压级 | 2 | 0.001 dB |
| 频带声压级 | 68 | 0.01 dB |
| 分声源频谱 | 476 | 0.01 dB |
| 分节点声级 | 180 | 0.05 dB |

同时检查时间轴、有限数值、步数及输出次数。容差涵盖不同求根／Newton 收敛处理和 Fortran 文本输出精度；不要求字节一致。`reference/results/*/reference.json` 保存源提交、输入及结果哈希、参考程序哈希。9 m/s 参考程序只增加诊断输出，物理计算仍全部为原 Fortran。

### 声学公式直接对照

```sh
cmake -S . -B build -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build -j 4
python tests/compare_fortran.py build/libaeroacoustics_shared.so --fortran-compiler gfortran --report build/kernels.json
```

`tests/fortran_reference.py` 从 `reference/fortran/` 的原声学例程建立独立 Fortran 动态库；原数值体不改公式，以双精度编译。测试包括 204 组参数、16 组 TNO 参数，共 56,576 个值。Python 负责生成参数、调用两个库及统计误差，没有 Python 声学数值实现。

### 重新计算原 Fortran 基线

这一步仅供验证开发，不是运行 C++ 算例的前置条件。需要 GFortran、CMake、BLAS/LAPACK。下载脚本使用固定 OpenFAST 与 r-test 提交；在新下载的原版源码目录执行以下命令。

Linux 示例：

```sh
python tools/fetch_openfast.py
cmake -S _work/openfast -B _work/build-reference -DCMAKE_BUILD_TYPE=Release -DDOUBLE_PRECISION=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build _work/build-reference --target openfast -j 4
python tools/run_full_case.py --executable _work/build-reference/glue-codes/openfast/openfast --label fortran-reference
python tests/compare_full_case.py _work/runs/fortran-reference build/official-check --report build/recomputed-reference.json
```

Windows 使用 MinGW，并按 [工具链说明](development.md#toolchain) 将 BLAS/LAPACK 指向 MKL。`run_full_case.py` 为原程序创建独立目录，保存日志、输入和程序哈希，确认正常结束；已有目录不会被覆盖。无需 C++ 接入补丁。

### 气动、结构和网格模块

先编译上一步的原 Fortran 库以及 C++ 测试程序。下面使用 Linux 默认 BLAS/LAPACK；Windows 在结构脚本末尾加 `--link-library 'MKL 的 mkl_rt.lib 路径'`，在 UA 脚本末尾加 `--mkl-lib '同一路径'`，并将 MKL DLL 目录加入 PATH。

```sh
python tests/compare_structural_modules.py --openfast-build _work/build-reference --cpp-build build --work _work/module-check --report build/modules.json
./build/turbine_module_probe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/ua-cpp.csv
python tests/compare_turbine_modules.py --openfast-build _work/build-reference --case examples/IEA_LB_RWT-AeroAcoustics --cpp-csv build/ua-cpp.csv --work _work/ua-check
```

| 对照范围 | 样本数 | 本机最大绝对差异 |
| --- | ---: | ---: |
| 声学公式，包括 TNO | 56,576 个数值 | 4.55e−13 dB |
| UnsteadyAero 升力／阻力／力矩系数 | 36,000 行 | 4.76e−14 |
| 结构加速度 | 5,700 行 | 7.96e−13 m/s² |
| 网格运动及载荷映射 | 4,700 行 | 7.67e−11，见各通道量纲 |
| BEM 残差、诱导系数与损失因子 | 7,200 行 | 6.31e−12 |

结构探针使用独立给定的模态位移、速度和载荷；网格探针给 Fortran 与 C++ 提供相同的运动和载荷；BEM 探针给定入流角，直接比较诱导公式。整机回归另覆盖求根、UA 状态传递及耦合时间推进。

若需追踪逐步状态，`tests/instrument_solver_reference.py` 可向开发用原源码加入 `solver-reference.csv` 诊断输出；重新编译后，配合 `turbine_solver_probe` 比较模态状态、气动状态、节点位置与载荷。该工具保存原文件副本，不改变求解公式。

```sh
python tests/compare_solver_trace.py _work/runs/fortran-trace/solver-reference.csv build/solver-cpp.csv --report build/solver-trace.json
```

逐步轨迹对照每个工况覆盖 288,000 行、6,624,000 个状态及载荷值。官方工况的入流角最大差异为 2.43e−6 rad，模态位移最大差异为 1.11e−5 m；最大局部升力系数差异约 9.90e−4，局部力差异约 3.053 N/m。后两项出现在非定常模型分段附近，未描述为逐位相同。详见 [官方轨迹](validation-solver-trace.json) 和 [9 m/s 轨迹](validation-solver-trace-wind9.json)。

### 运行依赖检查

Windows 的 `tests/check_runtime_windows.py` 递归检查 PE 导入表，将可执行文件、三个 C++ 运行库和官方输入复制到新目录，将子进程 PATH 限制为系统目录，然后运行完整算例并检查结果。

```powershell
python tests/check_runtime_windows.py build/aeroacoustics_turbine.exe _work/runtime-check --gcc-bin D:/Code_Configuration/gcc-16.2.0/mingw64/bin --report build/runtime.json
```

本机普通构建未导入 Fortran、OpenFAST、Python 或 MKL 动态库。结果见 [validation-runtime.json](validation-runtime.json)。这种检查证明已测试可执行文件的依赖情况，不等同于对任意后端或未来修改的保证。

### 图表与记录

```sh
python tools/plot_validation.py reference/results/IEA_LB_RWT-AeroAcoustics build/official-check
```

绘图另需 Matplotlib。曲线读取实际输出；平均频谱先转声能求平均，再转换为 dB。导出文件为 `examples/results/official_oaspl.csv`、`official_mean_spectrum.csv` 和 `docs/full-case-validation.png`。

当前整机报告为 [官方工况](validation-full-case.json) 与 [9 m/s 工况](validation-wind9.json)，模块报告为 [普通声学库](validation-fortran-portable.json)、[MKL 声学库](validation-fortran-mkl.json)、[UA](validation-unsteady-aero.json)、[结构／网格／BEM](validation-turbine-modules.json)。测试检验程序与原模型的一致程度，不包含实测噪声标定。

<a id="reliability"></a>

## 声级、文件输出与查表诊断（R1–R3）

本次处理非有限声级、零声能输出语义、文件写入检查和查表范围诊断。物理公式及默认端点插值策略保持不变。

### R1：区分静音与计算异常

- `-Inf` 继续表示零声能，参与总量求和时贡献为零。
- `NaN`、`+Inf`、声级转声能时溢出，以及累加得到的非有限/负声能均报错。
- 声源模型产生非法声级时，整机错误信息包含时间、观察点、叶片、原始节点编号、机制、频率和尾缘/入流模型开关。
- 单截面 C++ 接口报告机制与频率；工作区额外报告观察点和节点。独立库用户未指定诊断编号时，编号为 0，另提供工作区内的节点序号。

关键实现：[acoustic_levels.hpp](../include/acoustic_levels.hpp)、[spectrum.cpp](../src/acoustics/spectrum.cpp)、[workspace.cpp](../src/acoustics/workspace.cpp)、[driver.cpp](../src/acoustics/driver.cpp)。异常会使整机返回非零退出码，不再把它当作静音继续输出。

### R2：保持旧文件，同时给出零声能掩码

原有四类 `.out` 保留 Fortran 对照格式，其中零声能仍使用 `0 dB` 占位。每份实际生成的 `.out` 都增加同名 `.mask`，行、时间和通道顺序完全对应：

| mask | 含义 | 如何解释 `.out` |
| --- | --- | --- |
| 0 | 当前聚合结果为零声能 | `0 dB` 是占位值；可来自关闭的机制、未选中的节点或数值上为零的贡献 |
| 1 | 当前聚合结果为正声能 | 声级是有效结果；即使它恰好为 `0 dB`，也不是静音占位 |

掩码描述声能是否为零，不表示预测可信度，也不细分零声能的物理原因。非法值会报错，不会写成 mask=0。每份文件仍有三行表头；`run.json` 记录该编码规则。

[CheckedOutput](../include/checked_output.hpp) 对声学输出、掩码、动力学、查表报告和运行记录统一检查打开、写入、刷新与关闭。错误包含失败文件路径；完成提示只在所有输出关闭成功后出现。

开始生成新输出时先清空原 `run.json`，在数据文件均完成后才写入新运行记录。失败可能留下部分数据文件；调用者必须检查退出码，不能仅以 `.out` 存在判断成功。输入/参数解析阶段失败属于尚未开始输出的新运行，应同样以退出码为准。

### R3：默认汇总诊断，可选严格报错

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

### 验证与复现

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

<a id="interfaces"></a>

## 接口、状态与并行

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build -j 4
./build/acoustic_workspace_probe verify build/workspace.csv
./build/turbine_workspace_probe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst
./build/mechanism_probe
./build/simulation_api_probe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/api-check
./build/semantics_modules_probe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst
./build/performance_api_probe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/performance-api
python tests/check_semantics.py build/aeroacoustics_turbine build/semantics-cli
python tests/check_performance.py build build/performance-cli
```

- 声源测试覆盖 216 组开关组合、标准通道顺序和 TNO 替代关系。
- 状态测试覆盖模型所有权、复制/移动、检查点重放、重置、失败恢复，以及输入文件移走后已构造实例的重置。帧引用的有效期见 [仿真接口](development.md#simulation)。
- 量纲测试覆盖 PSD 积分、频带合并、声压/声功率换算、A 计权及非法频带；积分器使用含非对角耦合的 1、2、3、4、6 自由度系统核对解析响应。
- 性能接口测试检查完整/分块频谱、按需输出、scaled 求解器和独立工况并行。预热后 50 个默认结构步的 `operator new/new[]` 调用数为零；该统计不含初始化、文件输出、异常及新增工程模式。
- GitHub Actions 在 Linux 上执行数值与异常检查，并用 ThreadSanitizer 检查批量和工程模型的并发路径；不检测 MKL 动态库内部。

历史报告：[S1–S3](validation-structure.json)、[S4–S5](validation-semantics.json)、[P1–P3](validation-coupling.json)、[P4–P7](validation-performance.json)。报告各自记录其提交和测试环境，不将不同阶段的耗时拼成连续加速比。

<a id="engineering"></a>

## 非稳态控制与传播（E1–E2）

```sh
./build/engineering_probe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst examples/engineering/controller.dat
python tests/check_engineering.py build/aeroacoustics_turbine build/engineering-cli
```

C++ 探针检查四维风场插值、两质量传动链解析响应、变转速/变桨/偏航运动学导数、转矩映射守恒、空气吸收、刚性地面 6.0206 dB 极限及屏障组合。还检查稳态/无损退化、完整状态重放和共享只读配置的并发一致性。

CLI 检查正常与降噪模式、状态/执行器约束、四档主步长和风场越界拒绝。示例使用 `DT=0.0015625 s`；在 2 秒步长检查中，它相对 `0.00078125 s` 的最大总声级差为约 0.0211 dB，不能外推为所有控制参数及长时间工况的精度保证。

2026-09-19 的普通与 MKL 结果、程序哈希和默认工况回归见 [validation-engineering.json](validation-engineering.json)。新增模型尚无同配置 Fortran 或实测对照，控制器参数未作机型标定。模型方程、输入和适用范围见 [工程模型](engineering.md)。
