# 源码优化与性能验证

本次优化针对结构载荷重复映射、声源重复计算、临时数组分配、固定配置重复解析，以及经验模型中沿用的 Fortran 式组织方式。基线为 `ad44bfd`，物理模型、输入文件、时间网格和收敛阈值保持不变。

## 修改内容

| 问题 | 修改 | 主要代码 |
| --- | --- | --- |
| 数值 Jacobian 扰动一片叶片时，重新映射三片叶片的载荷 | 新增单叶片映射函数；每轮 Newton 迭代由 30 次单叶片映射减少为 12 次 | [rotor.cpp](../src/turbine/coupling/rotor.cpp)、[solver.cpp](../src/turbine/coupling/solver.cpp) |
| 每个观察点重复计算同一节点的边界层、声源谱形和 TNO 积分 | 分为 `prepare_*()` 与 `emit_*()`；先按节点准备声源，再按观察点施加距离和指向性 | [source_models.hpp](../src/acoustics/source_models.hpp)、[spectrum.cpp](../src/acoustics/spectrum.cpp)、[workspace.cpp](../src/acoustics/workspace.cpp) |
| 时间循环、频谱装配和 TNO 积分反复创建或复制数组 | 工作区保存节点、频谱和积分数组；主程序复用四类输出数组，并通过只读快照指针读取结果 | [driver.cpp](../src/acoustics/driver.cpp)、[tno.cpp](../src/acoustics/tno.cpp)、[turbine_main.cpp](../src/apps/turbine_main.cpp) |
| 节点和输出循环反复解析固定输入，重复校验参数和计算 A 计权 | 初始化时缓存 `SkewOptions`、`NrOutFile`、已校验声学参数和 A 计权；动态截面输入仍逐次检查 | [rotor.cpp](../src/turbine/coupling/rotor.cpp)、[workspace.cpp](../src/acoustics/workspace.cpp) |
| 经验公式集中在一个大文件，参数列表长，变量声明和数组下标沿用 Fortran 风格 | 按边界层、BPM 尾缘、其他 BPM 声源、入流、谱形和指向性拆分；使用具名数据类型、局部变量及从零开始的循环 | [源码索引](../src/README.md)、[kernel_compat.cpp](../src/acoustics/kernel_compat.cpp) |

BPM 的共享边界层在每个节点、每个声学采样时刻只计算一次。TNO 两侧积分也只准备一次，观察点变化只影响后续距离和指向性计算。缓存不会跨时间保留过期声源：新的速度、攻角、湍流强度和边界层输入会重新生成当前谱形。

经验公式的系数、分段条件、钝度频谱逐频带累积归一化和声级下限沿用原实现。此次没有增加噪声模型或改变气动与结构的推进次序。

## 性能结果

环境：Windows 11、Intel Core i7-11800H、GCC 16.2.0、Release（`-O3 -DNDEBUG`），普通 C++ 后端，未启用 MKL。基线和优化版先预热，再交替运行五轮，以下取中位数。原始样本、可执行文件 SHA256 和输出差异保存在 [validation-optimization.json](validation-optimization.json)。

### 完整 20 秒算例

| 工况 | 优化前 | 优化后 | 加速比 | 耗时减少 |
| --- | ---: | ---: | ---: | ---: |
| 官方 8 m/s | 3.536 s | 1.824 s | 1.94× | 48.4% |
| 9 m/s 扰动 | 3.497 s | 1.794 s | 1.95× | 48.7% |

计时来自 `run.json`，覆盖程序内的初始化、3200 步推进、201 次声学采样及输出处理，不包含外部启动进程的开销。没有缩短物理时长或减少输出。

### 声学局部测试

每次快照使用 30 个节点和 34 个频率；BPM 测量批次包含 20 次快照，TNO 包含 2 次，表内折算为一次快照。基线调用原 `snapshot_spectrum()`，优化版复用 `AcousticWorkspace`，因此加速同时包含声源共享和缓冲区复用的收益。

| 模型 | 观察点数 | 优化前 | 优化后 | 加速比 |
| --- | ---: | ---: | ---: | ---: |
| BPM | 2 | 2.013 ms | 0.659 ms | 3.06× |
| BPM | 16 | 16.508 ms | 1.169 ms | 14.13× |
| TNO | 2 | 739.289 ms | 361.894 ms | 2.04× |
| TNO | 16 | 6078.342 ms | 364.610 ms | 16.67× |

多观察点共享计算时收益更明显。以上是特定测试输入的局部结果，不能直接作为整机加速比。官方整机默认使用 BPM，不启用 TNO；此次数据也不是 MKL 与普通后端的速度比较。

## 数值检查

- **Fortran 公式对照**：普通和 MKL 后端分别通过 56,576 个值、204 组配置的对照，包含 16 组 TNO 配置。最大绝对差异约 `4.55×10⁻¹³ dB`，阈值为 `2×10⁻⁸ dB`。见 [普通后端报告](validation-fortran-portable.json) 和 [MKL 报告](validation-fortran-mkl.json)。
- **整机 Fortran 对照**：官方 8 m/s 和 9 m/s 扰动各检查 145,926 个声学值及 101 个输入文件，均满足原有阈值。见 [官方工况报告](validation-full-case.json) 和 [扰动工况报告](validation-wind9.json)。
- **优化前后整机对照**：两个工况的 `dynamics.csv` 均逐字节一致；9 m/s 的四份声学文件也逐字节一致。8 m/s 仅频带输出文件出现约 `1×10⁻¹⁰ dB` 的末位变化，其余三份一致。
- **工作区复用**：36 组配置与状态组合改变节点数、观察点数、攻角、TI、边界层和 A 计权，检查复用输出与独立截面调用，并与旧 C++ 库比较 21,420 个值；有限值最大差异约 `5.68×10⁻¹⁴ dB`。普通和 MKL 后端均通过复用检查。
- **声学驱动状态**：两个 TI 计算分支、边界层表、展向选点、厚度提取和积分辅助功能与旧 C++ 输出对照通过；驱动频谱最大差异约 `1.14×10⁻¹³ dB`。

公式对照误差与整机误差含义不同。整机仍包含 C++ 和原 Fortran 求根及耦合迭代的数值差异，原有支持范围见 [模块说明](standalone-port.md)。

## 接口与复现

`AcousticWorkspace::evaluate()` 返回内部快照的只读引用，`AcousticDriver::step_view()` 在采样时刻返回只读指针，其他时刻返回空指针。调用者须在下一次相应调用前使用完数据；需要持久保存时自行复制。各求解实例分别持有工作区，同一工作区不应被多个线程同时写入。

原 `AcousticDriver::step()` 继续返回独立快照；`section_spectrum()`、`snapshot_spectrum()`、标量模型函数和 C ABI 保留。应用使用新 C++ 头文件后应重新编译。整机仍可仅用 C++ 编译和运行；Python 脚本仅负责测试调度和数值文件比较。

构建与核心验证：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_BUILD_TESTS=ON
cmake --build build -j 4
./build/acoustic_workspace_probe verify build/workspace.csv
python tests/compare_fortran.py build/libaeroacoustics_shared.so --fortran-compiler gfortran --report build/kernels.json
python tests/run_standalone.py build/aeroacoustics_turbine build/official-check --report build/full-case.json
python tests/run_perturbation.py build/aeroacoustics_turbine build/wind9 --report build/wind9.json
```

在同一机器上构建优化前基线，并将新探针以兼容模式链接至旧库。下面使用单配置生成器、GCC/Clang：

```sh
git worktree add --detach build-baseline-src ad44bfd
cmake -S build-baseline-src -B build-baseline -DCMAKE_BUILD_TYPE=Release
cmake --build build-baseline -j 4
c++ -std=c++17 -O3 -DNDEBUG -DAERO_BENCHMARK_LEGACY -Ibuild-baseline-src/include tests/acoustic_workspace_probe.cpp build-baseline/libaeroacoustics.a -o build-baseline/acoustic_workspace_probe
python tests/benchmark_optimization.py build-baseline build --wind9-case build/wind9/inputs/IEA_LB_RWT-AeroAcoustics.fst --repeats 5 --baseline-commit ad44bfd --report build/optimization.json
```

Windows / MinGW 先运行 `. ./tools/use_gcc.ps1`，使用 `-G 'MinGW Makefiles'`，将 `.so` 改为 `libaeroacoustics_shared.dll` 并为可执行文件加 `.exe`；手动编译探针时将 `c++` 换为 `g++`。9 m/s 输入由前面的 `run_perturbation.py` 生成。复现实测值时应保持编译选项一致，并避免同时执行其他高负载任务。
