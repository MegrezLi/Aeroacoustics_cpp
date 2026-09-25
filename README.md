# Aeroacoustics_cpp

用 C++17 独立运行 OpenFAST 官方 `IEA_LB_RWT-AeroAcoustics` 整机声学算例。从 `.fst`、叶片、翼型和声学输入文件计算风场、叶片变形、气动力及噪声，完成 0–20 秒仿真。

默认编译和运行只需要 CMake 与 C++ 编译器。可执行程序不调用 OpenFAST、Fortran 或 Python，也不读取预先计算的节点时序。仓库内的 Fortran 源码和结果用于验证；Python 仅保留下载、测试、绘图及工具链安装脚本。

## 编译与运行

需要 CMake ≥3.20、支持 C++17 的编译器。官方算例的 101 个输入文件已经包含在仓库中。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j 4
./build/aeroacoustics_turbine examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/results
```

Windows / MinGW：

```powershell
. ./tools/use_gcc.ps1
cmake -S . -B build -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
./build/aeroacoustics_turbine.exe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/results
```

`use_gcc.ps1` 默认使用本机 GCC 16.2.0；其他位置可加 `-GccBin '你的 mingw64/bin 路径'`。使用 Visual Studio 多配置生成器时，可执行文件位于 `build/Release/`。已有不同生成器的构建目录应另取名称。

输出到指定目录：

| 文件 | 内容 |
| --- | --- |
| `IEA_LB_RWT-AeroAcoustics_1.out` | 两个观察点的总声压级 |
| `IEA_LB_RWT-AeroAcoustics_2.out` | 各观察点的 34 个频带 |
| `IEA_LB_RWT-AeroAcoustics_3.out` | 各频带的 7 类声源分量 |
| `IEA_LB_RWT-AeroAcoustics_4.out` | 各叶片、节点对观察点的声级贡献 |
| `dynamics.csv` | 每个时间步的 9 个模态位移和速度 |
| `run.json` | 运行统计、频带与计权、参考声压、模块配置及自由度布局 |
| 同名 `.mask` | 与每份 `.out` 对应：0 表示零声能占位，1 表示有效正声能 |
| `lookup_diagnostics.csv` | 按表格、坐标轴、叶片、节点和调用阶段统计查表超界 |

默认计算 3,200 个时间步，步长 0.00625 s；声学输出间隔 0.1 s，共 201 个时刻。声级参考声压为 20 μPa。程序另可接受第三个参数设置运行时长，例如末尾加 `2` 运行 2 秒。

非法声级和文件写入失败会使程序报错退出。`.out` 中的零声能仍使用原 `0 dB` 占位，须结合 `.mask` 与物理上的有效 `0 dB` 区分。查表默认保留端点值并汇总警告；可加 `--lookup-policy=error` 在首次超界时报错。字段、适用范围及测试说明见 [R1–R3 可靠性说明](docs/validation.md#reliability)。

## 模型范围

移植范围是该算例实际启用的物理模块：

- **InflowWind**：恒定风、幂律风切变及风向变换。
- **ElastoDyn**：三叶片，每片两个挥舞模态和一个摆振模态；包括结构扭转、质量、刚度、阻尼、重力、离心力、科里奥利项及轴向缩短。
- **AeroDyn / BEM**：叶素动量法、Prandtl 根尖损失、切向诱导、Buhl 高诱导修正、Pitt/Peters 偏斜尾流修正。
- **AirfoilInfo / UnsteadyAero**：翼型表线性插值、`UA_Mod=3` 的 Minnema/Pierce 非定常气动状态。
- **网格传递与耦合**：结构到气动节点的运动插值、气动到结构的载荷映射、广义 α 时间积分。
- **AeroAcoustics**：BPM 尾缘噪声、钝度、层流、叶尖、Lowson 入流噪声、Simplified Guidati 厚度修正，以及 TNO 模型和 61 点积分。

官方工况使用 8 m/s 风速、10.04 rpm 固定转速和 1.17° 固定桨距。新增工程配置可启用独立的两质量传动链、发电机转矩响应和闭环变桨/偏航；塔架与平台仍固定。这些新增模型不是 ServoDyn/ROSCO 的完整移植。

数值求根与结构 Newton 迭代采用 C++ 实现，收敛处理与原求解器有差别，不承诺逐位一致。具体支持项、求解顺序及输入限制见 [模块说明](docs/development.md#scope)。默认整机算例未启用 TNO，TNO 由单独的 Fortran 对照测试覆盖。

## 非稳态运行与传播

```sh
./build/aeroacoustics_turbine examples/engineering/closed_loop.fst build/engineering --wind-grid=examples/engineering/gust_veer.wind --controller=examples/engineering/controller.dat --propagation=examples/engineering/propagation.dat
```

风场接口读取时间/x/y/z 网格中的三分量速度；闭环控制器计算转速、桨距、偏航及正常/降噪模式响应。传播层增加均匀空气吸收、平面地面反射和主导山脊/薄屏障衍射，保留源模型已有的距离与指向性。

示例控制参数尚未针对真实机组标定，风场为确定性测试数据。新增模型已有解析和时间步收敛检查，尚无同配置 Fortran 或实测对照；传播模型不含气象梯度折射，也不等同于完整 ISO 9613-2。输入格式、参数单位、调用顺序与验证范围见 [E1–E2 工程说明](docs/engineering.md)。

## 表面状态与接收点统计

```sh
./build/aeroacoustics_turbine examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/metrics --surfaces=examples/engineering/surfaces.dat --metrics=examples/engineering/metrics.dat
```

`--surfaces` 接入有来源、有效攻角/Re 范围及不确定度说明的翼型/边界层数据；粗糙度、侵蚀和转捩标签不自动生成声级修正。`--metrics` 输出接收时间历程、LAeq、时间百分位、风速分箱、AM 描述指标和受声点地图，并可加入外部指定的运动音调。

示例数据仅用于测试。新增统计不等同于 IEC/IOA 完整测量流程；接收时延暂不与相干地面反射组合，窄带音调也不由现有宽带模型自动预测。参数与输出见 [工程模型](docs/engineering.md#surface-data)。

## 多机尾流与风场声学

```sh
./build/aeroacoustics_farm examples/farm/farm.dat build/farm-results
```

`aeroacoustics_farm` 按来流顺序独立求解各机组，用气动推力更新准稳态 Jensen 尾流，并在共同受声点汇总叶片声与外部给定的机械/冷却声功率。每台风机可以分别配置闭环控制、表面数据及固定塔架势流/塔影；单机也可用 `--tower=examples/farm/tower-test.dat` 启用塔架入流影响。

随附双机布局、塔架和机械声源均为测试数据。当前限于固定水平风向、对齐机组、平坦场址和远尾流，支持自由场/空气吸收；按共同源时间汇总，不含动态尾流、附加湍流和输运时延。输入格式、推力限幅诊断与适用范围见 [E6 工程说明](docs/engineering.md#farm)。

## Fortran 与 C++ 对比

参考源码固定为 OpenFAST [`2895884`](https://github.com/OpenFAST/openfast/tree/2895884d2be01862173c88d70f86b358d2f1a50a)，输入固定为 r-test [`dd5feaa`](https://github.com/OpenFAST/r-test/tree/dd5feaaaa500ba7283140107806300d551cff0a7/glue-codes/openfast/IEA_LB_RWT-AeroAcoustics)。原 Fortran 以双精度编译。

| 20 秒整机输出：最大绝对差异 | 官方 8 m/s | 9 m/s 扰动工况 |
| --- | ---: | ---: |
| 总声压级 | 0.000078 dB | 0.000016 dB |
| 频带声压级 | 0.000206 dB | 0.000082 dB |
| 分声源频谱 | 0.000852 dB | 0.000177 dB |
| 分节点声级 | 0.000767 dB | 0.005637 dB |

每个工况比较全部 145,926 个声学值，并核对 101 个输入文件的 SHA256。改变风速后，总声压级相对原工况最多改变约 1.526 dB。

模块对照另覆盖 56,576 个声学数值、36,000 组非定常翼型状态、5,700 行结构状态、4,700 行网格映射和 7,200 组 BEM 诱导系数。详细误差见 [验证与复现](docs/validation.md#fortran) 和 `docs/validation-*.json`。Windows 运行依赖检查只装入 C++ 运行库，在 PATH 仅含系统目录时完成了整机算例，见 [运行检查报告](docs/validation-runtime.json)。

![Fortran 与独立 C++ 整机对比](docs/full-case-validation.png)

以下测试需要 Python、NumPy；声学源码直接对比还需要 GFortran，它们不参与风机求解。

```sh
cmake -S . -B build -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build --config Release -j 4
python tests/compare_fortran.py build/libaeroacoustics_shared.so --fortran-compiler gfortran --report build/kernels.json
python tests/run_standalone.py build/aeroacoustics_turbine build/official-check --report build/full-case.json
python tests/run_perturbation.py build/aeroacoustics_turbine build/wind9-check --report build/wind9.json
```

Windows 将 `.so` 换成 `libaeroacoustics_shared.dll`，并给可执行文件加 `.exe`。GitHub Actions 执行 Linux 编译、这三项对比、声学与整机缓冲区检查，以及声源组合、库接口、检查点、工程控制/传播、并发和异常检查。已有 Fortran 源码与结果保存在 `reference/`。

## C++ 仿真接口

整机调度已独立为 `Simulation`。可以逐帧调用 `next()`，或通过 `run_case()` 运行完整工况；自定义 `ResultSink` 可直接接收内存中的结果。模型数据共享且只读，求解状态通过检查点完整保存和恢复。

声源配置集中管理了七类声源和模型替代关系。命令行及标准输出格式保持不变；旧接口中直接读取 `solver.state`、`solver.time` 的 C++ 调用需改为 `state()`、`time()` 等只读访问器。示例、检查点范围和迁移表见 [仿真接口说明](docs/development.md#simulation)。

声学接口明确区分 PSD、频带均方声压、声压级与声功率，保存频带上下限和计权状态；整机聚合拒绝将任意密集频率样本作为独立频带累加。A 计权输出标题为 dBA。通用广义 α 积分器按具名自由度与耦合块求解，当前三叶片后端通过加速度接口接入。接口、频带约定和扩展边界见 [S4–S5 说明](docs/development.md#quantities-modules)。

## 批量运行与求解选项

独立工况可通过 C++ 线程并行运行，输出目录必须互不重叠：

```sh
./build/aeroacoustics_batch 2 case8/IEA_LB_RWT-AeroAcoustics.fst results/wind8 case9/IEA_LB_RWT-AeroAcoustics.fst results/wind9
```

整机默认按一个观察点一块计算，并按 `NrOutFile` 分配和聚合输出。`--observer-block-size=N` 可调整块大小。`--solver=scaled` 启用尺度化扰动、残差判停和同一步内 Jacobian 复用；默认 `reference` 保留原数值路径。参数、限制和验证见 [P4–P7 说明](docs/performance.md#p4-p7)。

## 性能

时间循环复用结构、气动与声学工作区；节点声源计算在观察点间共享，输出按需聚合。已有优化的测试条件、计时与复现命令统一放在 [性能文档](docs/performance.md)。独立工况可并行运行；MKL 是否加速取决于 TNO 积分规模，不能据此推断整机一定更快。

## 截面示例与 MKL

`aeroacoustics_example` 保留单截面频谱示例：

```sh
./build/aeroacoustics_example spectrum.csv
```

Intel oneMKL 是可选后端，用于 TNO 积分的向量指数和 BLAS 求和。启用时配置 `-DAEROACOUSTICS_USE_MKL=ON` 和 `MKL_DIR`，运行时提供 MKL 动态库。MinGW 使用 `mkl_rt` 单一动态接口。安装及命令见 [工具链说明](docs/development.md#toolchain)。

## 独立验证与不确定性

`aeroacoustics_validate` 使用 C++ 计算测量对照误差、参数扰动敏感性、相关不确定性预算和等权样本分位区间；可直接导入整机 `receiver_map.csv` 的 LAeq。校准与验证批次分开统计，未知不确定度留空，不自动判定实机精度或标准符合性。

```sh
./build/aeroacoustics_validate compare examples/validation/observations.csv examples/validation/predictions.csv build/data-comparison 2
```

随附数据为合成测试数据。格式、完整命令、适用假设和真实数据接入步骤见 [E8 验证说明](docs/validation.md#independent-validation)。

## 文档导航

| 文档 | 内容 |
| --- | --- |
| [工作流](aeroacoustics工作流.md) | 输入参数、模块调用顺序和主要函数 |
| [开发接口](docs/development.md) | 源码索引、支持范围、C++ API、状态与工具链 |
| [工程模型](docs/engineering.md) | 非稳态风、闭环控制、室外传播、表面数据与多机尾流 |
| [验证](docs/validation.md) | Fortran 对照、解析检查、异常测试和复现命令 |
| [性能](docs/performance.md) | 分块、并行、求解选项和历史性能数据 |
| [待优化](待优化.md) | 已完成项、原始代码审查和后续工程需求 |

## 源码

```text
src/
├── acoustics/       # 经验声学模型、TNO、频谱装配、观察点几何和时间驱动
├── numerics/        # Gauss–Kronrod 积分、普通/MKL 数值后端
├── propagation/     # 空气吸收、地面镜像和屏障衍射
├── turbine/
│   ├── aerodynamics/ # BEM 与非定常气动
│   ├── inflow/       # 网格风场、塔架入流与准稳态多机尾流
│   ├── control/      # 两质量传动链、发电机和闭环执行器
│   ├── analysis/     # 接收时间、窄带线声源输运和工程统计
│   ├── structure/    # 叶片模态动力学
│   ├── coupling/     # 网格映射、转子装配和时间积分
│   ├── simulation/   # 仿真调度、声学适配、结果聚合与输出
│   └── io/           # 风机算例、翼型及风场输入
├── interfaces/      # C ABI
└── apps/            # 截面、整机、批量、风场与验证程序入口
```

输入参数与模块调用顺序见 [aeroacoustics 工作流](aeroacoustics工作流.md)，目录职责及源文件索引见 [源码索引](docs/development.md#source)。公开头文件仍位于 `include/`。

- `include/turbine/batch.hpp`：独立工况并行与逐项运行结果。
- `include/turbine/integrator.hpp`、`include/turbine/modules.hpp`：通用积分器、自由度布局和模块组合校验。
- `include/acoustic_quantities.hpp`：频带、声学量类型与单位换算。
- `include/turbine/simulation.hpp`：整机库接口、逐步运行与进程内检查点。
- `include/mechanisms.hpp`：声源机制、标准通道顺序及单位。
- `include/aeroacoustics.hpp`、`include/aeroacoustics_c.h`：声学库 C++ 接口与 C ABI。
- `examples/IEA_LB_RWT-AeroAcoustics/`：可直接运行的官方输入。
- `tests/`、`reference/`：Fortran 对照、回归测试与参考结果。
- `tools/`：下载、绘图及工具链辅助脚本。

Apache-2.0；来源及修改说明见 [NOTICE](NOTICE) 和 [LICENSE](LICENSE)。
