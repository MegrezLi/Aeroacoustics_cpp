# Fortran 与独立 C++ 验证

所有数值对照均以原 Fortran 为参考。默认 C++ 求解器不使用 `reference/`；测试脚本在运行结束后读取参考结果并比较。

## 随仓库运行回归

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

## 声学公式直接对照

```sh
cmake -S . -B build -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build -j 4
python tests/compare_fortran.py build/libaeroacoustics_shared.so --fortran-compiler gfortran --report build/kernels.json
```

`tests/fortran_reference.py` 从 `reference/fortran/` 的原声学例程建立独立 Fortran 动态库；原数值体不改公式，以双精度编译。测试包括 204 组参数、16 组 TNO 参数，共 56,576 个值。Python 负责生成参数、调用两个库及统计误差，没有 Python 声学数值实现。

## 重新计算原 Fortran 基线

这一步仅供验证开发，不是运行 C++ 算例的前置条件。需要 GFortran、CMake、BLAS/LAPACK。下载脚本使用固定 OpenFAST 与 r-test 提交；在新下载的原版源码目录执行以下命令。

Linux 示例：

```sh
python tools/fetch_openfast.py
cmake -S _work/openfast -B _work/build-reference -DCMAKE_BUILD_TYPE=Release -DDOUBLE_PRECISION=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build _work/build-reference --target openfast -j 4
python tools/run_full_case.py --executable _work/build-reference/glue-codes/openfast/openfast --label fortran-reference
python tests/compare_full_case.py _work/runs/fortran-reference build/official-check --report build/recomputed-reference.json
```

Windows 使用 MinGW，并按 [工具链说明](toolchain.md) 将 BLAS/LAPACK 指向 MKL。`run_full_case.py` 为原程序创建独立目录，保存日志、输入和程序哈希，确认正常结束；已有目录不会被覆盖。无需 C++ 接入补丁。

## 气动、结构和网格模块

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

## 运行依赖检查

Windows 的 `tests/check_runtime_windows.py` 递归检查 PE 导入表，将可执行文件、三个 C++ 运行库和官方输入复制到新目录，将子进程 PATH 限制为系统目录，然后运行完整算例并检查结果。

```powershell
python tests/check_runtime_windows.py build/aeroacoustics_turbine.exe _work/runtime-check --gcc-bin D:/Code_Configuration/gcc-16.2.0/mingw64/bin --report build/runtime.json
```

本机普通构建未导入 Fortran、OpenFAST、Python 或 MKL 动态库。结果见 [validation-runtime.json](validation-runtime.json)。这种检查证明已测试可执行文件的依赖情况，不等同于对任意后端或未来修改的保证。

## 图表与记录

```sh
python tools/plot_validation.py reference/results/IEA_LB_RWT-AeroAcoustics build/official-check
```

绘图另需 Matplotlib。曲线读取实际输出；平均频谱先转声能求平均，再转换为 dB。导出文件为 `examples/results/official_oaspl.csv`、`official_mean_spectrum.csv` 和 `docs/full-case-validation.png`。

当前整机报告为 [官方工况](validation-full-case.json) 与 [9 m/s 工况](validation-wind9.json)，模块报告为 [普通声学库](validation-fortran-portable.json)、[MKL 声学库](validation-fortran-mkl.json)、[UA](validation-unsteady-aero.json)、[结构／网格／BEM](validation-turbine-modules.json)。测试检验程序与原模型的一致程度，不包含实测噪声标定。
