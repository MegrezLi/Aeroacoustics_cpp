# Aeroacoustics_cpp

C++17 气动噪声程序：读取 OpenFAST 格式输入，计算叶片结构、气动力和各观察点的噪声频谱。仓库提供可直接运行的 20 秒整机算例，以及复现 NASA BPM 1989 表 D2 的单翼段算例。

编译和运行需要 **CMake ≥3.20、C++17 编译器**。默认程序独立运行；Fortran 用于参考对照，Python 用于测试和绘图。

## 编译并运行整机算例

所有命令在仓库根目录执行。官方算例的 101 个输入文件已包含在仓库中。

Windows / MinGW，在 PowerShell 中执行：

```powershell
. ./tools/use_gcc.ps1
cmake -S . -B build -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
./build/aeroacoustics_turbine.exe examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/results
```

`use_gcc.ps1` 默认查找 `D:/Code_Configuration/gcc-16.2.0/mingw64/bin`。编译器装在其他位置时，在第一行加 `-GccBin '你的 mingw64/bin 路径'`。

Linux / macOS，使用单配置生成器：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
./build/aeroacoustics_turbine examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/results
```

Visual Studio 构建时加 `--config Release`，可执行文件位于 `build/Release/`。更换生成器时使用新的构建目录。

默认风速 **8 m/s**、转速 **10.04 rpm**、桨距 **1.17°**；运行 20 s，时间步长 0.00625 s，声学间隔 0.1 s，共 201 个声学时刻。只试跑 2 s 时，在运行命令末尾加 `2`。

## 查看结果

结果写入 `build/results`，声压级以 **20 μPa** 为参考，默认未加 A 计权。

| 文件 | 内容 |
| --- | --- |
| `IEA_LB_RWT-AeroAcoustics_1.out` | 两个观察点的总声压级 |
| `IEA_LB_RWT-AeroAcoustics_2.out` | 各观察点的 34 个频带 |
| `IEA_LB_RWT-AeroAcoustics_3.out` | 各频带的 7 类声源分量 |
| `IEA_LB_RWT-AeroAcoustics_4.out` | 各叶片、节点对观察点的声级贡献 |
| 同名 `.mask` | 0：零声能占位；1：有效正声能，与 `.out` 逐项对应 |
| `dynamics.csv` | 每步的 9 个模态位移和速度 |
| `run.json` | 运行统计、声学单位、频带、计权及模块配置 |
| `lookup_diagnostics.csv` | 查表超界的位置、次数和输入范围 |

读取 `.out` 时同时读取 `.mask`：静音占位和有效声级都可能写成 `0 dB`。计算或写入失败返回非零退出码。查表默认使用端点值并记录越界；运行命令加 `--lookup-policy=error` 可在首次越界时停止。详见[输出与诊断](docs/validation.md#reliability)。

## 选择算例和功能

先学声学公式，从 **[BPM 单翼段教程](examples/bpm1989/README.md)** 开始；查看整机数据流，读 **[整机工作流](aeroacoustics工作流.md)**。

| 任务 | 入口与配置 | 说明 |
| --- | --- | --- |
| 复现 BPM 1989 | `aeroacoustics_bpm1989` | NACA 0012，27 个频带，与原程序表 D2 对照 |
| 变风速、变桨、变转速 | `aeroacoustics_turbine` + `--wind-grid`、`--controller` | [网格风场、两质量传动链和闭环控制](docs/engineering.md#wind-control) |
| 空气吸收、地面、屏障 | `--propagation` | [均匀介质直线路径传播](docs/engineering.md#propagation) |
| 更换翼型或边界层 | `--surfaces` | [数据格式、两侧顺序、攻角和 Re 范围](docs/engineering.md#surface-data) |
| LAeq、百分位、AM、测点地图 | `--metrics` | [接收时间统计及外部运动音调](docs/engineering.md#receiver-metrics) |
| BPM / TNO / Howe–Chase 对比 | `--trailing-edge` | [直尾缘、锯齿及 TNO 外缘速度模式](docs/engineering.md#trailing-edge) |
| 多机尾流和声能汇总 | `aeroacoustics_farm` | [Jensen 尾流、固定塔架、机械声源输入](docs/engineering.md#farm) |
| 多个独立工况并行 | `aeroacoustics_batch` | [工况并行、观察点分块及求解选项](docs/performance.md) |
| 预测与测量配对、不确定性 | `aeroacoustics_validate` | [compare、budget、ensemble 与数据格式](docs/validation.md#independent-validation) |

工程示例中的风场、控制参数、表面数据和机械声谱用于软件测试；应用到具体机组时替换为对应输入。表 D2 和原 OpenFAST 结果用于模型实现对照，测量数据通过 `aeroacoustics_validate` 单独比较。

## 模型与接口

整机参考路径包括稳态风、三叶片九模态结构、BEM、`UA_Mod=3` 的 Minnema/Pierce 非定常气动、运动/载荷网格映射及广义 α 积分。塔架与平台固定；闭环控制使用项目独立的控制器。支持的输入组合见[模块范围](docs/development.md#scope)。

声学固定输出七个通道：`LBL`、`TBL_pressure`、`TBL_suction`、`TBL_separation`、`bluntness`、`tip`、`inflow`。TNO 替换压力面和吸力面尾缘声，保留 BPM 分离声；Howe–Chase 替换两侧尾缘声并关闭分离声。各模型通过[机制配置](docs/development.md#simulation)选择，避免重复叠加。

嵌入 C++ 程序时用 `turbine::run_case()`；逐步处理结果用 `Simulation::next()`，自定义输出用 `ResultSink`。状态复制、检查点、声学量及源码索引见[开发接口](docs/development.md)。MKL 声学后端的构建方法见[工具链](docs/development.md#toolchain)。

## 运行测试

需要 Python 3.11+ 和 NumPy。使用上面已配置的 `build` 目录：

```sh
python -m pip install 'numpy>=1.26,<3'
cmake -S . -B build -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build --config Release -j 4
ctest --test-dir build -C Release --output-on-failure -j 2
```

测试包括单元、命令行、BPM 算例及已保存的 8/9 m/s 整机 Fortran 对照。要现场编译 Fortran 声学参考库，再配置 `-DAEROACOUSTICS_TEST_FORTRAN=ON`，并确保 GFortran 在 PATH 中。

整机对照每个工况比较 145,926 个声学值。已保存记录中的最大绝对差值为：

| 输出 | 8 m/s | 9 m/s |
| --- | ---: | ---: |
| 总声级 | 0.000078 dB | 0.000016 dB |
| 频带 | 0.000206 dB | 0.000082 dB |
| 分声源 | 0.000852 dB | 0.000177 dB |
| 分节点 | 0.000767 dB | 0.005637 dB |

容差、参考版本、模块对照和重算方法见[验证与复现](docs/validation.md)。原始报告：[8 m/s](docs/validation-full-case.json)、[9 m/s](docs/validation-wind9.json)。

## 文档与源码

| 文档 | 内容 |
| --- | --- |
| [BPM 教程](examples/bpm1989/README.md) | 编译、文献参数、公式、调用模块和参数练习 |
| [整机工作流](aeroacoustics工作流.md) | 输入引用、时间推进、BEM/UA、结构和声学调用链 |
| [工程模型](docs/engineering.md) | 风场、控制、传播、表面数据、统计、尾流及尾缘模型 |
| [开发接口](docs/development.md) | 源码索引、C++ API、状态与单位、工具链 |
| [验证与复现](docs/validation.md) | CTest、Fortran 对照、测量配对与不确定性 |
| [性能与求解选项](docs/performance.md) | 缓存、分块、并行、Newton 配置与计时 |

公开接口在 `include/`，实现按声学、气动、结构、耦合和应用分布于 `src/`。测试在 `tests/`，运行输入在 `examples/`，固定 Fortran 参考源码与结果在 `reference/`。

授权与来源见 [LICENSE](LICENSE) 和 [NOTICE](NOTICE)。
