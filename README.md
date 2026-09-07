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
| `run.json` | 步长、运行时长、步数和声学输出次数 |

默认计算 3,200 个时间步，步长 0.00625 s；声学输出间隔 0.1 s，共 201 个时刻。声级参考声压为 20 μPa。程序另可接受第三个参数设置运行时长，例如末尾加 `2` 运行 2 秒。

## 模型范围

移植范围是该算例实际启用的物理模块：

- **InflowWind**：恒定风、幂律风切变及风向变换。
- **ElastoDyn**：三叶片，每片两个挥舞模态和一个摆振模态；包括结构扭转、质量、刚度、阻尼、重力、离心力、科里奥利项及轴向缩短。
- **AeroDyn / BEM**：叶素动量法、Prandtl 根尖损失、切向诱导、Buhl 高诱导修正、Pitt/Peters 偏斜尾流修正。
- **AirfoilInfo / UnsteadyAero**：翼型表线性插值、`UA_Mod=3` 的 Minnema/Pierce 非定常气动状态。
- **网格传递与耦合**：结构到气动节点的运动插值、气动到结构的载荷映射、广义 α 时间积分。
- **AeroAcoustics**：BPM 尾缘噪声、钝度、层流、叶尖、Lowson 入流噪声、Simplified Guidati 厚度修正，以及 TNO 模型和 61 点积分。

官方工况使用 8 m/s 风速、10.04 rpm 固定转速和 1.17° 固定桨距。塔架、传动链、平台及控制器在此算例中关闭；这些自由度和模块不属于本项目的移植范围。程序会拒绝未支持的主要模型开关，不能直接替代 OpenFAST 运行任意风机。

数值求根与结构 Newton 迭代采用 C++ 实现，收敛处理与原求解器有差别，不承诺逐位一致。具体支持项、求解顺序及输入限制见 [模块说明](docs/standalone-port.md)。默认整机算例未启用 TNO，TNO 由单独的 Fortran 对照测试覆盖。

## Fortran 与 C++ 对比

参考源码固定为 OpenFAST [`2895884`](https://github.com/OpenFAST/openfast/tree/2895884d2be01862173c88d70f86b358d2f1a50a)，输入固定为 r-test [`dd5feaa`](https://github.com/OpenFAST/r-test/tree/dd5feaaaa500ba7283140107806300d551cff0a7/glue-codes/openfast/IEA_LB_RWT-AeroAcoustics)。原 Fortran 以双精度编译。

| 20 秒整机输出：最大绝对差异 | 官方 8 m/s | 9 m/s 扰动工况 |
| --- | ---: | ---: |
| 总声压级 | 0.000078 dB | 0.000016 dB |
| 频带声压级 | 0.000206 dB | 0.000082 dB |
| 分声源频谱 | 0.000852 dB | 0.000177 dB |
| 分节点声级 | 0.000767 dB | 0.005637 dB |

每个工况比较全部 145,926 个声学值，并核对 101 个输入文件的 SHA256。改变风速后，总声压级相对原工况最多改变约 1.526 dB。

模块对照另覆盖 56,576 个声学数值、36,000 组非定常翼型状态、5,700 行结构状态、4,700 行网格映射和 7,200 组 BEM 诱导系数。详细误差见 [验证与复现](docs/full-case.md) 和 `docs/validation-*.json`。Windows 运行依赖检查只装入 C++ 运行库，在 PATH 仅含系统目录时完成了整机算例，见 [运行检查报告](docs/validation-runtime.json)。

![Fortran 与独立 C++ 整机对比](docs/full-case-validation.png)

以下测试需要 Python、NumPy；声学源码直接对比还需要 GFortran。它们不参与风机求解。

```sh
cmake -S . -B build -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build --config Release -j 4
python tests/compare_fortran.py build/libaeroacoustics_shared.so --fortran-compiler gfortran --report build/kernels.json
python tests/run_standalone.py build/aeroacoustics_turbine build/official-check --report build/full-case.json
python tests/run_perturbation.py build/aeroacoustics_turbine build/wind9-check --report build/wind9.json
```

Windows 将 `.so` 换成 `libaeroacoustics_shared.dll`，并给可执行文件加 `.exe`。GitHub Actions 执行 Linux 编译和这三项对比。已有 Fortran 源码与结果保存在 `reference/`，旧 Python 算法及 Python 数值基线已移除。

## 截面示例与 MKL

`aeroacoustics_example` 保留单截面频谱示例：

```sh
./build/aeroacoustics_example spectrum.csv
```

Intel oneMKL 是可选后端，用于 TNO 积分的向量指数和 BLAS 求和。启用时配置 `-DAEROACOUSTICS_USE_MKL=ON` 和 `MKL_DIR`，运行时提供 MKL 动态库。MinGW 使用 `mkl_rt` 单一动态接口。安装及命令见 [工具链说明](docs/toolchain.md)。

## 源码

- `src/turbine/`、`include/turbine/`：输入、结构、BEM、非定常气动、网格映射及整机求解器。
- `src/kernels.cpp`、`src/aeroacoustics.cpp`、`src/driver.cpp`：声学公式、积分和时间状态。
- `include/aeroacoustics.hpp`、`include/aeroacoustics_c.h`：声学库 C++ 接口与 C ABI。
- `examples/IEA_LB_RWT-AeroAcoustics/`：可直接运行的官方输入。
- `tests/`、`reference/`：Fortran 对照、回归测试与参考结果。
- `tools/`：下载、绘图及工具链辅助脚本。

Apache-2.0；来源及修改说明见 [NOTICE](NOTICE) 和 [LICENSE](LICENSE)。这些测试验证代码移植，不构成风机实测噪声验证。
