# 验证与复现

[Fortran 对照](#fortran) · [声级、输出与查表](#reliability) · [接口与并行](#interfaces) · [工程模型](#engineering) · [独立数据与不确定性](#independent-validation)

先运行 CTest 检查程序，再按需求重算 Fortran 基线或导入测量数据。命令均在仓库根目录执行；手动调用测试脚本时使用新的输出目录。Fortran 对照衡量移植一致性，测量对照衡量实际预测误差。

## 统一测试入口

需要 Python 3.11+ 与 NumPy；直接编译 Fortran 对照还需 GFortran。

```sh
python -m pip install 'numpy>=1.26,<3'
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_BUILD_TESTS=ON -DAEROACOUSTICS_TEST_FORTRAN=ON
cmake --build build -j 4
ctest --test-dir build --output-on-failure -j 2
```

Windows/MinGW 先执行 `. ./tools/use_gcc.ps1`，配置加 `-G "MinGW Makefiles"`。Visual Studio 等多配置生成器：构建加 `--config Release`，CTest 加 `-C Release`。没有 GFortran 时省略 `AEROACOUSTICS_TEST_FORTRAN=ON`，仍可对照已保存的整机 Fortran 结果。

| 标签 | 内容 |
| --- | --- |
| `unit` | 声源开关、C ABI、声能合成、工作区、状态、统计、模型与输入检查 |
| `integration` | 各命令行及组合模型回归 |
| `reference` | 8/9 m/s 整机基线；启用时增加声学内核和尾缘 Fortran 对照 |
| `example` | 截面示例、BPM 文献算例及分析脚本 |
| `concurrency` | 五组含并发调用的原生测试，Linux CI 另用 ThreadSanitizer 执行 |

例如 `ctest --test-dir build -L unit --output-on-failure` 只运行单元检查。每次执行的输出保留在 `build/test-results/测试名-随机后缀/`，可重复运行且不会互相覆盖。失败日志会给出具体目录。完整原 OpenFAST 的结构/UA 对照需要另行构建上游程序，命令见下文。

声源测试同时覆盖 NaN/+Inf 拒绝、−Inf 静音、C ABI 失败时输出不变；统计测试覆盖部分静音百分位、移动声源能量守恒及不同 ReceiverDT 下的积分一致性。

<a id="fortran"></a>

## Fortran 与独立 C++ 对照

默认 C++ 求解器不读取 `reference/`。回归脚本在运行后读取 Fortran 结果比较；测量数据通过[独立数据验证](#independent-validation)导入。

### 随仓库运行回归

先按 README 编译。以下命令在项目根目录执行，需要 Python 3.11+ 与 NumPy。Windows 对可执行文件加 `.exe`。

```sh
python tests/check_reference.py build/aeroacoustics_turbine build/official-check --report build/full-case.json
python tests/check_reference.py --wind-speed 9 build/aeroacoustics_turbine build/wind9-check --report build/wind9.json
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

配置时启用 `AEROACOUSTICS_TEST_FORTRAN=ON`，然后运行：

```sh
ctest --test-dir build -C Release -R "^fortran_kernels$" --output-on-failure
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
python tests/compare_full_case.py _work/runs/fortran-reference build/official-check/outputs --report build/recomputed-reference.json
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

逐步轨迹对照每个工况覆盖 288,000 行、6,624,000 个状态及载荷值。官方工况的入流角最大差异为 2.43e−6 rad，模态位移最大差异为 1.11e−5 m；最大局部升力系数差异约 9.90e−4，局部力差异约 3.053 N/m。后两项出现在非定常模型分段附近。详见 [官方轨迹](validation-solver-trace.json) 和 [9 m/s 轨迹](validation-solver-trace-wind9.json)。

### 运行依赖检查

Windows 的 `tests/check_runtime_windows.py` 递归检查 PE 导入表，将可执行文件、三个 C++ 运行库和官方输入复制到新目录，将子进程 PATH 限制为系统目录，然后运行完整算例并检查结果。

```powershell
python tests/check_runtime_windows.py build/aeroacoustics_turbine.exe _work/runtime-check --gcc-bin D:/Code_Configuration/gcc-16.2.0/mingw64/bin --report build/runtime.json
```

归档普通构建仅依赖系统及 C++ 运行库，未导入 Fortran、OpenFAST、Python 或 MKL 动态库。程序哈希与依赖结果见 [validation-runtime.json](validation-runtime.json)。更换后端后重新运行上述检查。

### 图表与记录

```sh
python tools/plot_validation.py reference/results/IEA_LB_RWT-AeroAcoustics build/official-check/outputs
```

绘图另需 Matplotlib。曲线读取实际输出；平均频谱先转声能求平均，再转换为 dB。导出文件为 `examples/results/official_oaspl.csv`、`official_mean_spectrum.csv` 和 `docs/full-case-validation.png`。

当前整机报告为 [官方工况](validation-full-case.json) 与 [9 m/s 工况](validation-wind9.json)，模块报告为 [普通声学库](validation-fortran-portable.json)、[MKL 声学库](validation-fortran-mkl.json)、[UA](validation-unsteady-aero.json)、[结构／网格／BEM](validation-turbine-modules.json)。

<a id="reliability"></a>

## 声级、文件输出与查表诊断

以下说明非有限声级、零声能输出语义、文件写入检查和查表范围诊断。物理公式及默认端点插值策略保持不变。

### 静音与计算异常

公开 `db_sum()` 在求最大值前检查每个声级，拒绝 NaN/+Inf；C ABI 按操作校验频率和所用物理字段，并在复制任何通道前检查全部输出。失败时返回非零、保留调用方缓冲区，不增加成功计数；未使用字段与 TNO 的有符号外缘速度比/非正 Cf 关闭语义保持兼容。`mechanism_probe` 覆盖声能合成与 C ABI，`metrics_probe` 覆盖部分静音百分位、不同接收网格、移动声源、裁剪窗口、背景声能和风速分箱的统计守恒。两者均通过 CTest 接入 CI。

- `-Inf` 继续表示零声能，参与总量求和时贡献为零。
- `NaN`、`+Inf`、声级转声能时溢出，以及累加得到的非有限/负声能均报错。
- 声源模型产生非法声级时，整机错误信息包含时间、观察点、叶片、原始节点编号、机制、频率和尾缘/入流模型开关。
- 单截面 C++ 接口报告机制与频率；工作区额外报告观察点和节点。独立库用户未指定诊断编号时，编号为 0，另提供工作区内的节点序号。

关键实现：[acoustic_levels.hpp](../include/acoustic_levels.hpp)、[spectrum.cpp](../src/acoustics/spectrum.cpp)、[workspace.cpp](../src/acoustics/workspace.cpp)、[driver.cpp](../src/acoustics/driver.cpp)。异常会使整机返回非零退出码，不再把它当作静音继续输出。

### 零声能掩码与文件写入

原有四类 `.out` 保留 Fortran 对照格式，其中零声能仍使用 `0 dB` 占位。每份实际生成的 `.out` 都增加同名 `.mask`，行、时间和通道顺序完全对应：

| mask | 含义 | 如何解释 `.out` |
| --- | --- | --- |
| 0 | 当前聚合结果为零声能 | `0 dB` 是占位值；可来自关闭的机制、未选中的节点或数值上为零的贡献 |
| 1 | 当前聚合结果为正声能 | 声级是有效结果；即使它恰好为 `0 dB`，也不是静音占位 |

掩码描述声能是否为零，不表示预测可信度，也不细分零声能的物理原因。非法值会报错，不会写成 mask=0。每份文件仍有三行表头；`run.json` 记录该编码规则。

[CheckedOutput](../include/checked_output.hpp) 对声学输出、掩码、动力学、查表报告和运行记录统一检查打开、写入、刷新与关闭。错误包含失败文件路径；完成提示只在所有输出关闭成功后出现。

开始生成新输出时先清空原 `run.json`，在数据文件均完成后才写入新运行记录。失败可能留下部分数据文件；调用者必须检查退出码，不能仅以 `.out` 存在判断成功。输入/参数解析阶段失败属于尚未开始输出的新运行，应同样以退出码为准。

### 查表诊断与严格模式

默认策略 `clamp`：超出表格范围后仍使用端点值，并在成功运行后生成 `lookup_diagnostics.csv`；有越界时终端给出一次汇总警告。

```sh
./build/aeroacoustics_turbine examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/results --lookup-policy=clamp
```

严格策略 `error`：首次超界即失败，错误信息包含文件、坐标轴、输入值、表格范围、时间、叶片、节点和查询阶段。

```sh
./build/aeroacoustics_turbine examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/strict 20 --lookup-policy=error
```

时长参数与策略参数可以交换顺序；Windows 可执行文件加 `.exe`。严格模式也检查 BEM 试探角度和非定常翼型内部查询，因此严格模式失败并不自动说明最终工作攻角越界。

报告按表格、坐标轴、叶片、节点和阶段分组，包含：

- `calls`：实际超界查询次数；同一步重复查询也计数，一个 BL 查询的两个轴越界分别计数。
- `first_time_s`、`last_time_s`：该组超界查询的最早与最晚输入时刻，不是持续越界时长。
- `min_value`、`max_value`、`lower`、`upper`：超界输入极值和表格有效范围。
- `axis`：`alpha_deg`（度）或 `Re`（无量纲）；翼型攻角先沿用原来的周期归一化，再检查表格范围。
- `stage`：`BEM_trial`、`UA_evaluate`、`UA_advance` 或 `BL_interpolate`。它们用于区分求根试探和不同调用阶段，不能直接当作最终物理状态的超界比例。

精确落在端点属于有效输入，不计越界；没有越界时仍生成仅含表头的报告。失败时使用终端定位信息，部分或空的 CSV 不能作为完整统计。

库级入口保持原签名。需要报告时，可使用 [lookup_diagnostics.hpp](../include/lookup_diagnostics.hpp) 的 `LookupReport`、`LookupSession` 和 `LookupLocation`。会话是显式启用、线程局部的，嵌套结束后恢复调用者上下文；未启用会话的库调用保持端点策略。多线程用户应为每个线程建立独立报告，再合并结果。

### 验证与复现

```sh
ctest --test-dir build -C Release -R "^(reliability|cli_reliability|reference_wind[89])$" --output-on-failure
```

CTest 自动创建独立工作目录；手动运行 CLI 脚本时指定新目录。测试仅修改目录中的输入副本。

针对性检查覆盖有效 0 dB 与静音、NaN/+Inf、声能溢出、模型异常上下文、真实查表的端点/超界、严格模式、嵌套诊断会话、输出路径冲突及刷新失败。Linux 还使用 `/dev/full` 检查真实写入失败。

归档结果见 [validation-reliability.json](validation-reliability.json)：普通/MKL 可靠性探针、工作区和各 56,576 个声学值的 Fortran 对照通过；8/9 m/s 的 `.out` 与 `dynamics.csv` 与 `3c24ded` 逐字节一致。CLI 覆盖输出路径冲突、非有限声级及严格查表。

<a id="interfaces"></a>

## 接口、状态与并行

```sh
ctest --test-dir build -C Release -R "^(acoustic_workspace|turbine_workspace|mechanisms|simulation_api|semantics|performance_api|cli_semantics|cli_performance)$" --output-on-failure
```

- 声源测试覆盖 216 组开关组合、标准通道顺序和 TNO 替代关系。
- 状态测试覆盖模型所有权、复制/移动、检查点重放、重置、失败恢复，以及输入文件移走后已构造实例的重置。帧引用的有效期见 [仿真接口](development.md#simulation)。
- 量纲测试覆盖 PSD 积分、频带合并、声压/声功率换算、A 计权及非法频带；积分器使用含非对角耦合的 1、2、3、4、6 自由度系统核对解析响应。
- 性能接口测试检查完整/分块频谱、按需输出、scaled 求解器和独立工况并行。预热后 50 个默认结构步的 `operator new/new[]` 调用数为零；该统计不含初始化、文件输出、异常及新增工程模式。
- GitHub Actions 在 Linux 上执行数值与异常检查，并用 ThreadSanitizer 检查批量和工程模型的并发路径；不检测 MKL 动态库内部。

历史报告：[仿真与状态](validation-structure.json)、[量纲与模块](validation-semantics.json)、[耦合缓存](validation-coupling.json)、[分块与并行](validation-performance.json)。各报告记录对应提交和测试环境。

<a id="engineering"></a>

## 非稳态控制与传播

```sh
ctest --test-dir build -C Release -R "^(engineering|cli_engineering)$" --output-on-failure
```

C++ 探针检查四维风场插值、两质量传动链解析响应、变转速/变桨/偏航运动学导数、转矩映射守恒、空气吸收、刚性地面 6.0206 dB 极限及屏障组合。还检查稳态/无损退化、完整状态重放和共享只读配置的并发一致性。

CLI 检查正常/降噪模式、执行器约束、四档主步长及风场越界。固定 2 s 工况中，`DT=0.0015625 s` 相对 `0.00078125 s` 的最大总声级差约 0.0211 dB。

普通/MKL 结果与程序哈希见 [validation-engineering.json](validation-engineering.json)。输入方程与组合限制见[工程模型](engineering.md)。

<a id="surface-metrics"></a>

## 表面数据、接收时间与统计

```sh
ctest --test-dir build -C Release -R "^(metrics|cli_metrics)$" --output-on-failure
```

C++ 解析测试检查非均匀时间积分、持续时间百分位、恒定声能、静音、正弦调制、匀速接近与静止音调、跨风速箱积分、边界层量纲/两侧顺序、有效范围和表观声功率几何关系。状态测试检查工程历史重放、重置、共享只读表面输入下的并行运行及内存上限。

整机 CLI 检查：仅替换边界层时动力学逐字节保持一致而声学发生变化；提供原完整翼型文件时重现同一状态；未计权与已 A 计权输入得到相同 LAeq；减少普通输出数量不改变统计；风速箱时间与能量之和回到总体统计；独立音调产生可解析的接收频率变化。另检查四点地图、背景声能相加、时间窗和不兼容配置拒绝，以及闭环风场/空气吸收的联合运行。

固定测试工况中，`DT_AA/ReceiverDT` 从 0.1 s 缩小到 0.05 s，四点 LAeq 最大差约 0.00108 dB。输入、区间和程序哈希见 [validation-metrics.json](validation-metrics.json)。

<a id="independent-validation"></a>

## 独立数据验证、敏感性与不确定性

`aeroacoustics_validate` 提供测量配对、局部不确定性预算和等权样本统计。随附 CSV 为合成测试数据；评估实际预测误差时，用风洞或现场数据替换测量表，并填写来源与工况。工具不自动拟合或修正模型。

### 运行示例

以下目录必须尚不存在，父目录须已存在；Windows 为程序名添加 `.exe`。

```sh
./build/aeroacoustics_validate compare examples/validation/observations.csv examples/validation/predictions.csv build/data-comparison 2
./build/aeroacoustics_validate budget examples/validation/perturbations.csv examples/validation/correlations.csv build/uncertainty-budget 2
./build/aeroacoustics_validate ensemble examples/validation/samples.csv build/ensemble-summary
```

比较和预算命令的末尾参数是覆盖因子 k；`k=2` 的覆盖率取决于分布和自由度。三个命令保存输入副本，全部输出成功后写入 `report.json`。

### 测量与预测的配对

两份 CSV 使用相同列名，列顺序可变，字段如下。格式示例见 [observations.csv](../examples/validation/observations.csv) 和 [predictions.csv](../examples/validation/predictions.csv)。

| 列 | 定义 |
| --- | --- |
| `id`、`receiver` | 配对记录唯一编号、接收点编号 |
| `group`、`split` | 独立试验/采集批次标识；`calibration` 或 `validation`。同一批次的相关时段必须使用同一 group，不能跨两种划分 |
| `origin`、`provenance` | 测量侧为 `field`、`wind_tunnel` 或 `synthetic`；预测侧为 `simulation`。来源说明记录数据文件、仪器/处理方法或计算提交、输入与运行结果位置 |
| `quantity`、`weighting` | `LAeq` 必须为 A；`band_Leq` 可为 A 或 Z。声压基准均为 20 μPa，不接受声功率或 PSD 混入 |
| `frequency_hz`、`bandwidth_hz` | 频带中心和带宽；LAeq 两者填写 0，频带声级两者必须为正。两侧须采用相同频带定义、滤波方法与积分频率范围，并在来源中注明 |
| `signal` | `turbine_only` 或 `total`，分别表示风机声或包含背景的总声；本程序不代做实测背景扣除 |
| `start_s`、`end_s` | 对应的平均时窗。先在外部统一时间原点；窗口必须相同，不自动插值、补零或外推 |
| `x_m`、`y_m`、`z_m` | 同一坐标系下的接收点位置 |
| `wind_mps`、`direction_deg`、`surface` | 一致定义的代表风速、[0,360) 风向与表面状态标识；标准化风速、气象代表性和工况筛选由数据准备过程确定 |
| `level_db`、`standard_u_db` | 声级和 dB 单位的标准不确定度。未知不确定度留空，不得用 0 代替未知；0 表示该项确实按零处理 |

配对必须覆盖两侧全部 ID，无重复或缺失。上下文字符串须一致，数值容差为 `1e-8*max(1,abs(a),abs(b))`，仅用于浮点文本差异。无效数据须在前处理阶段筛选，并在来源中保留筛选依据。

`residuals.csv` 给出预测减测量的 dB 误差。两侧标准不确定度已知时，按模型误差与测量误差独立的假设计算 `u_diff=hypot(u_model,u_measurement)`。仅当 `u_diff>0` 时输出归一化残差及 `abs(error)<=k*u_diff` 标记；未知时相关列留空，不影响原始误差统计。

`groups.csv` 按划分、来源、试验批次、量、计权、频带、信号、表面状态、1 m/s 风速箱和 30° 风向扇区分组，输出记录数、偏差、MAE、RMSE、最大绝对误差及归一化统计。每条记录等权，在 dB 域统计；不据此计算均值置信区间或自动修正校准偏差。输入方负责声明数据来源和采集批次。

### 直接导入整机接收点结果

```sh
./build/aeroacoustics_validate import-map build/run/receiver_map.csv context.csv predictions.csv "C++ commit and input/run provenance"
```

`context.csv` 使用上述预测格式，`origin=simulation`，`level_db` 可先填 0。此命令按 `receiver` 读取 `receiver_map.csv`，核对测点坐标和实际时窗，按 `signal` 选风机 LAeq 或含背景 LAeq，并写入新的预测 CSV。仅支持 A 计权 LAeq；`standard_u_db` 保留上下文中的值。风况、表面状态、频率积分范围及运行版本需由使用者对照实际算例填写，导入器不会从测量值推断这些条件。频带测量目前使用统一 CSV 接口，尚无自动频带历史导入器。

### 参数敏感性和相关不确定性

[perturbations.csv](../examples/validation/perturbations.csv) 每行指定一个 `output_id` 下的一项参数：`parameter,unit,x_minus,x0,x_plus,u_x,y_minus_db,y0_db,y_plus_db,provenance`。输出可以是某接收点 LAeq 或某频带声级，但同一 output_id 必须代表同一输出条件。先通过整机/批量接口运行基准与逐参数正负扰动工况；其余参数、测点和统计设置保持一致，再填入真实计算结果。扰动步长由调用方设置。

要求正负扰动对称且所有参数使用同一基准输出，最多 256 项。计算中心差分 `c_i=(y_plus-y_minus)/(x_plus-x_minus)`，另报告左右单边导数、曲率及 `abs(c_i)*u_i`。建议至少用两种扰动步长、充分长的统计时窗核对稳定性；左右导数和曲率差异大时，局部线性传播可能不适用。

[correlations.csv](../examples/validation/correlations.csv) 列为 `output_id,parameter_a,parameter_b,rho`。每对参数只列一次，遗漏的对明确按不相关处理；空预算也须保留表头。相关系数必须在 [-1,1]，矩阵必须半正定（分解数值容差 1e-12），允许完全相关/反相关的奇异矩阵，拒绝不一致矩阵。相关性必须有数据或建模依据。

按 `u_y² = sum(c_i*c_j*u_i*u_j*rho_ij)` 计算标准不确定度，输出 `k*u_y` 和基准值上下界。各单项贡献不是存在相关性时可直接相加的方差份额。`u_x` 是各自参数单位下的标准不确定度，不是扰动步长、上下限或扩展不确定度。本实现采用一阶传播，关系式依据 [NIST 不确定性传播说明](https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-appendix-law-propagation-uncertainty)。输出区间只包含预算中列出的误差项。

### 等权样本结果

[samples.csv](../examples/validation/samples.csv) 列为 `output_id,sample_id,level_db,provenance`。输入外部按已声明联合分布抽样、经过 C++ 求解得到的等权结果；采样分布、参数相关性、随机种子、模型/输入版本和失败工况处理应记录在来源或其指向的记录中。本工具只汇总外部等权样本；规则参数扫描不能直接作为概率样本。

同一输出内 sample_id 唯一，各输出须拥有完全相同的样本集合，至少两个样本。`ensemble.csv` 输出 dB 均值、样本标准差、线性插值的 2.5/50/97.5 百分位及稳定声能平均值。分位区间描述输入样本的离散程度，不是均值置信区间；扩大样本量后检查分位数是否稳定。

### 验证命令与当前证据

```sh
ctest --test-dir build -C Release -R "^(validation|cli_validation)$" --output-on-failure
```

C++ 解析检查包含已知偏差/RMSE、未知或零不确定度、相关/反相关与非法矩阵、二次函数曲率和已知样本分位数。CLI 检查划分泄漏、量与频带不匹配、重复 ID、未知不确定度、输出复用和缺失样本的拒绝。

整机接入测试运行 7.8、7.9、8、8.1、8.2 m/s 五个 2 s 工况，核对接收结果导入和两种扰动步长的敏感度。记录见 [validation-e8.json](validation-e8.json)。

CSV 支持 UTF-8、可选 BOM、CRLF/LF、引号内逗号和双引号转义；不接受多行字段、空行、未知/重复列、非有限声级。每表最多 100,000 条记录，每行最多 1 MiB。静音或检出限截断数据需单独处理，目前不以任意低声级代替。

<a id="farm"></a>

## 塔架、尾流与多机声能汇总

```sh
ctest --test-dir build -C Release -R "^(farm|cli_farm)$" --output-on-failure
```

C++ 探针核对圆柱势流、Powles 塔影、Jensen 解析值、历史插值、方向和平移、气动力积分、点声源距离/指向性、音调吸收及状态/并发一致性。Python 生成输入并比较 C++ 输出。

CLI 检查单机退化（普通声学输出逐字节一致，LAeq 与原接收点统计一致）、下游风速与动力学反馈、上游结果不变、反转机组输入顺序后的汇总一致性、整体平移和 A 计权一致性。另检查各声源 LAeq 按声能相加回到风场总量，以及双机闭环控制、表面数据、塔架、空气吸收和固定冷却声源的联合运行。内存上限、过近尾流布局、重名和不支持的地面组合均须拒绝。

归档记录包括普通/MKL 测试、8/9 m/s 整机及声学 Fortran 对照、双机 20 s 完整输出，见 [validation-farm.json](validation-farm.json)。Linux CI 运行风场回归和 ThreadSanitizer 并发检查。

2 s 联合测试记录 72 个 CT 限幅样本，20 s 双机示例记录 73 个。限幅仅作用于 Jensen 尾流输入，不改变气动力。风场、塔架和机械谱为合成输入；模型条件见[多机说明](engineering.md#farm)。

<a id="trailing-edge"></a>

## 尾缘模型与外缘速度模式

```sh
ctest --test-dir build -C Release -R "^(trailing_edge|cli_trailing_edge|fortran_trailing_edge)$" --output-on-failure
```

先启用 `AEROACOUSTICS_BUILD_TESTS=ON`；其中 `fortran_trailing_edge` 还要求 `AEROACOUSTICS_TEST_FORTRAN=ON`。原 `TBLTE_TNO` 函数体保持不变，适配层传入非单位速度比。`tests/howe_reference.f90` 独立实现论文公式与复步长求导，用于数值对照。

测试覆盖：h=0 和趋零极限、直尾缘谱形、高频锯齿渐近、PSD/频带解析积分、距离及展长缩放、指向性零点、无机制重复计数；TNO reference/input、单位/非单位/带符号输入、旧内核兼容、非法输入，以及标量/工作区、检查点/重置和并发一致性。

整机检查包括相同边界层下更换声学模型不改变动力学、其他机制不变、观察点分块一致，以及闭环控制、塔架、表面数据、空气吸收、统计和风场入口组合。普通/MKL 构建和原 8/9 m/s 的 20 秒 Fortran/C++ 回归分别执行。Howe 直尾缘与锯齿另运行 20 秒合成整机对比。

结果与程序哈希见 [validation-trailing-edge.json](validation-trailing-edge.json)，输入及适用条件见[尾缘模型](engineering.md#trailing-edge)。Howe 的 20 s 对比使用合成边界层；实际降噪量需与相应测量数据配对。
