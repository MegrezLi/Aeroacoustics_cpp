# 从一个翼段开始：复现 BPM 1989 自噪声算例

本教程用一个固定 NACA 0012 翼段复现 BPM 1989 表 D2：编译 C++ 程序，生成 27 个频带的声压级，再核对文献结果。输入只有翼段、来流和观察点参数。

## 先把程序跑起来

在仓库根目录执行。需要 CMake 3.20 或更高版本，以及支持 C++17 的编译器。这个目标不依赖 Fortran、MKL 或 Python。

Linux / macOS，使用单配置生成器：

```sh
cmake -S . -B build-bpm -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_BUILD_EXAMPLES=ON -DAEROACOUSTICS_BUILD_TURBINE=OFF
cmake --build build-bpm --target aeroacoustics_bpm1989 -j 4
./build-bpm/aeroacoustics_bpm1989 build-bpm/comparison.csv --check
```

Windows / MinGW，在 PowerShell 中执行：

```powershell
# 如果 g++、mingw32-make 已在 PATH 中，可以省略这一行。
# 否则请把路径换成自己电脑上的 mingw64/bin。
. ./tools/use_gcc.ps1 -GccBin 'D:/Code_Configuration/gcc-16.2.0/mingw64/bin'
cmake -S . -B build-bpm -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_BUILD_EXAMPLES=ON -DAEROACOUSTICS_BUILD_TURBINE=OFF
cmake --build build-bpm --target aeroacoustics_bpm1989 -j 4
./build-bpm/aeroacoustics_bpm1989.exe build-bpm/comparison.csv --check
```

使用 Visual Studio 生成器时，构建命令加 `--config Release`，运行路径改为 `build-bpm/Release/aeroacoustics_bpm1989.exe`。不要在同一个构建目录内切换生成器。

`--check` 要求表 D2 的最大绝对差值不超过 **0.01 dB**，超出返回状态码 2。输出路径的父目录须已存在；不传参数时，在当前目录写入 `bpm1989_comparison.csv`。

仓库附带的结果见 [results/comparison.csv](results/comparison.csv)。2026-10-08 使用 GCC 16.2.0、Release 配置运行得到：

```text
Max absolute difference: 0.00620027 dB
RMSE: 0.00559309 dB
```

## 文献与对照数据

算例来自 Brooks、Pope 和 Marcolini 的经典报告：

> Brooks, T. F., Pope, D. S., and Marcolini, M. A. (1989). *Airfoil Self-Noise and Prediction*. NASA Reference Publication 1218.

[NASA 文献记录](https://ntrs.nasa.gov/citations/19890016302) · [报告全文](https://ntrs.nasa.gov/api/citations/19890016302/downloads/19890016302.pdf)

BPM 是三位作者姓氏的缩写。选取 Figure 45(a) 对应的 Appendix D 算例，因为报告同时列出了输入和逐频带输出。

| 要查什么 | 报告印刷页码 | PDF 阅读器页码（从 1 开始） |
| --- | ---: | ---: |
| Appendix D 示例输入 | 113 | 121 |
| Table D2：本算例的计算输出（27 个总声级已逐行对照原表图像） | 114 | 122 |
| 原程序的声速、运动黏度默认值 | 116 | 124 |

**表 D2 是原程序的预测值，本例比较这组模型输出。** Figure 45(a) 另含实验数据，不在本次对照中。

## 模型、声源与输出单位

空气流过翼型时，壁面附近会形成边界层。边界层内的扰动经过尾缘，能够向外辐射声音。BPM 用实验归纳的边界层关系和谱形函数来估计这种声音，不会在这里求解瞬时涡结构或声波传播方程。

本例打开了两组模型。湍流边界层尾缘噪声（TBL-TE）给出压力面、吸力面和分离相关的三个分量；层流边界层涡脱落噪声（LBL-VS）给出另一个分量。分离相关分量是 BPM 分解中的一项，它在小迎角下有数值，并不等于代码判定整个翼型已经失速。

`itrip=0` 表示不强制转捩，允许计算层流涡脱落分量。本例按原文同时开启 `ITURB` 和 `ILAM`，并非假定整个边界层始终为层流。

观察点处的结果是 **1/3 倍频程频带声压级**，参考声压为 20 μPa，没有施加 A 计权。中心频率沿用原表的 27 个名义频率，从 100 Hz 到 40000 Hz。不要把这些值当成 dB/Hz 的功率谱密度，也不要插入大量密集频点后直接累加为总声级。

## 设置文献参数

输入集中在 [main.cpp](main.cpp) 的 `Parameters p` 和 `Section s` 中。前者规定声学模型和空气参数，后者描述这个翼段与观察点。

| 物理量 | 文献设置 | 代码字段 |
| --- | ---: | --- |
| 翼型 | NACA 0012 | 使用 BPM 经验边界层关系，无坐标文件输入 |
| 弦长 c | 0.3048 m | `s.chord` |
| 展长 L | 0.4572 m | `s.span` |
| 来流速度 U | 71.3 m/s | `s.speed` |
| 有效迎角 α* | 1.516° | `s.alpha_deg` |
| 到观察点的距离 r | 1.22 m | `s.trailing.distance` |
| 指向角 θ、φ | 均为 90° | `s.trailing.theta`、`s.trailing.phi` |
| 声速 c₀ | 340.46 m/s | `p.spdsound` |
| 运动黏度 ν | 1.4529 × 10⁻⁵ m²/s | `p.kinvisc` |

迎角用附录中的有效值 `ALPSTAR=1.516°`；图中风洞几何迎角不能直接替代。声速和运动黏度按上表显式赋值。

对应的开关是：

```cpp
p.x_blmethod = 1;     // 根据弦长、Re、迎角和转捩条件估计边界层
p.itrip = 0;          // 不强制转捩
p.tbltemod = 1;       // BPM 湍流边界层尾缘模型
p.lammod = 1;         // 层流涡脱落模型
p.timod = 0;          // 不叠加入流湍流噪声
p.tipmod = 0;         // 不叠加翼尖涡噪声
p.bluntmod = 0;       // 不叠加钝尾缘噪声
p.aweighting = false;
```

`s.is_tip=false` 表示普通翼段，`s.tabulated_boundary_layer=false` 使用经验边界层。`s.stall_deg=12.5°` 是库的失速分支阈值，本例的 1.516° 迎角低于该值。

由这些参数可以得到两个基本量：

```math
M=\frac{U}{c_{0}}\approx0.2094,\qquad Re_{c}=\frac{Uc}{\nu}\approx1.496\times10^6.
```

前者进入声源尺度和指向性修正，后者用于边界层及谱形的经验关系。

## 声压级是怎样一步步算出来的

先估计边界层。`thick()` 根据弦长、雷诺数、迎角和转捩条件，计算压力面的边界层厚度 $\delta_{p}$，以及压力面、吸力面的位移厚度 $\delta_{p}^{\ast}$、$\delta_{s}^{\ast}$。层流模型使用前者；本例的湍流尾缘模型使用后两者。边界层厚度与位移厚度不是同一个量，接外部数据时尤其不能互换。

接着把频率写成无量纲的 Strouhal 数。尾缘模型使用

```math
St_{p}=\frac{f\delta_{p}^{\ast}}{U},\qquad St_{s}=\frac{f\delta_{s}^{\ast}}{U},
```

层流模型则使用压力面的边界层厚度 $\delta_{p}$ 构造频率尺度。这样的处理把不同速度、边界层厚度对应的谱形联系起来：特征长度或速度改变，谱峰的位置也会变化。

在当前未进入大分离分支的条件下，尾缘噪声的计算结构可概括为：

```math
\begin{aligned}
L_{p} &= 10\log_{10}\left(\frac{\delta_{p}^{\ast}M^5LD_{h}}{r^2}\right)+A_{p}+K_{1}-3+\Delta K_{1},\\
L_{s} &= 10\log_{10}\left(\frac{\delta_{s}^{\ast}M^5LD_{h}}{r^2}\right)+A_{s}+K_{1}-3,\\
L_{\alpha} &= 10\log_{10}\left(\frac{\delta_{s}^{\ast}M^5LD_{h}}{r^2}\right)+B+K_{2}.
\end{aligned}
```

这里的 A、B 描述谱形，K 系列是经验幅值修正；下标 p、s 分别表示压力面和吸力面。表达式沿用 BPM 的 SI 单位经验标定，不能只留下对数项就把它当作完整的声压级公式。原始分段系数较多，阅读时直接对照 [bpm_trailing_edge.cpp](../../src/acoustics/bpm_trailing_edge.cpp) 更容易确认每一项的来源。

层流涡脱落分量采用类似结构：

```math
L_{LBL}=10\log_{10}\left(\frac{\delta_{p} M^5LD_{h}}{r^2}\right)+G_{1}+G_{2}+G_{3}.
```

G₁ 控制相对峰值频率的谱形，G₂ 随雷诺数比变化，G₃ 包含迎角相关的幅值项。这些关系在 [bpm_other_sources.cpp](../../src/acoustics/bpm_other_sources.cpp) 的 `prepare_laminar()` 中实现。

上述 Dₕ 是朝向观察点的辐射修正。当前调用的是 `directh_te()`，其关系为

```math
D_{h}=\frac{2\sin^2(\theta/2)\sin^2\phi}{(1+M\cos\theta)[1+(M-M_{c})\cos\theta]^2},\qquad M_{c}=0.8M.
```

本例 θ=φ=90°，Dₕ 接近 1。代码为保持与参考实现一致，角度转换系数使用 `0.017453`；它不是精确的 π/180。距离衰减已经包含在声源公式的 r² 中，因此这里不用再额外叠加一次球面扩散。

最后，把各分量的声能相加：

```math
L_{total}(f)=10\log_{10}\left[\sum_{j}10^{L_{j}(f)/10}\right].
```

对应函数是 `db_sum()`。它先减去最大声级再求和，以改善数值稳定性。不要把四个 dB 数值直接相加；也要区分“同一频带内合并声源”与“跨频带求总体声级”。本例 CSV 的 `total_dB` 指前者。

## 调用模块

程序入口只需要调用一次：

```cpp
const auto spectra = section_spectrum(p, s);
```

实际经过的函数如下。`prepare` 阶段准备与观察点无关的边界层、谱形和幅值系数；`emit` 阶段加入观察点的距离与指向性。

```mermaid
flowchart TD
    Main[main.cpp: Parameters + Section] --> API[section_spectrum]
    API --> Select[SourceSelection / validate]
    API --> Prep[prepare_section]
    Prep --> BL[boundary_layer / thick]
    BL --> TE[prepare_trailing_edge]
    BL --> LBL[prepare_laminar]
    TE --> Shape[amin / amax / bmin / bmax / a0comp]
    API --> Emit[emit_section]
    Emit --> ET[emit_trailing_edge]
    Emit --> EL[emit_laminar]
    ET --> Dir[directh_te]
    EL --> Dir
    ET --> Channels[Mechanisms: 7 source spectra]
    EL --> Channels
    Channels --> Sum[main.cpp: db_sum per band]
    Sum --> CSV[CSV and Table D2 comparison]
```

| 文件 | 本例中的职责 |
| --- | --- |
| [main.cpp](main.cpp) | 设置工况，调用声学库，汇总频带声能，与表 D2 比较并写 CSV |
| [aeroacoustics.hpp](../../include/aeroacoustics.hpp) | 定义 `Parameters`、`Section`、`Geometry` 以及公开函数 |
| [spectrum.cpp](../../src/acoustics/spectrum.cpp) | 校验输入、选择模型，组织准备和输出阶段，提供 `db_sum()` |
| [boundary_layer.cpp](../../src/acoustics/boundary_layer.cpp) | `boundary_layer()`、`thick()`：经验边界层厚度 |
| [bpm_trailing_edge.cpp](../../src/acoustics/bpm_trailing_edge.cpp) | 准备并输出压力面、吸力面和分离相关的尾缘分量 |
| [bpm_spectral_shapes.cpp](../../src/acoustics/bpm_spectral_shapes.cpp) | 尾缘模型使用的分段谱形函数 |
| [bpm_other_sources.cpp](../../src/acoustics/bpm_other_sources.cpp) | 本例仅调用其中的层流涡脱落准备与输出函数 |
| [directivity.cpp](../../src/acoustics/directivity.cpp) | `directh_te()`：尾缘声源的观察方向修正 |
| [mechanisms.hpp](../../include/mechanisms.hpp) | 定义七类声源的枚举、输出次序和名称 |

本例执行 BPM 尾缘和层流涡脱落模型。`TnoWorkspace` 按需分配，`tbltemod=1` 时不分配 TNO 积分数组；整机气动、结构和时间推进模块不参与此算例。

## 读取频谱与误差

CSV 每行对应一个频带。`LBL_dB` 是层流分量；`TBL_pressure_dB`、`TBL_suction_dB`、`TBL_separation_dB` 是三个尾缘分量。后面的 `total_dB`、`reference_D2_dB`、`difference_dB` 分别给出当前计算、文献模型值和二者之差。

`bluntness_dB`、`tip_dB`、`inflow_dB` 都是 `-inf`，表示禁用声源的零声能。它不是程序溢出。原表这些列中的 `0.000` 是占位符，不能当成真实的 0 dB 参与叠加。另一种常见数值是尾缘模型保留的 `-100 dB` 下限；它是极小的有限声级，与禁用的 `-inf` 含义不同。

归档结果中的部分频带：

| 频率 / Hz | 当前计算 / dB | 表 D2 / dB | 差值 / dB |
| ---: | ---: | ---: | ---: |
| 1000 | 53.562617 | 53.568 | -0.005383 |
| 3150 | 60.490411 | 60.496 | -0.005589 |
| 4000 | 59.449186 | 59.455 | -0.005814 |
| 16000 | 49.584800 | 49.591 | -0.006200 |

![总声压级和差值](results/comparison.png)

两组总谱的峰值均在 **3150 Hz**，最大绝对差值在 **16000 Hz**。当前计算整体比表 D2 低约 0.0056 dB，偏差来源尚未定位。下图单独显示逐频带差值。

这里的 RMSE 是逐频带 dB 差值的均方根：

```math
\mathrm{RMSE}=\sqrt{\frac{1}{27}\sum_{i=1}^{27}(L_{i}-L_{D2,i})^2}.
```

它描述两组谱值的数值差异，不是宽带声压误差，也不是测量不确定度。

## 重新生成图和核查报告

`analyze.py` 会检查频率列表、有效声源的有限性、禁用声源的零能量、声能叠加以及差值。只生成 JSON 报告时，Python 标准库就够了：

```sh
python examples/bpm1989/analyze.py build-bpm/comparison.csv --output-dir build-bpm/report
```

绘图还需要 Matplotlib，建议在自己的 Python 虚拟环境中安装：

```sh
python -m pip install matplotlib
python examples/bpm1989/analyze.py build-bpm/comparison.csv --output-dir build-bpm/report --plot
```

加 `--plot` 时输出 `summary.json`、`comparison.png` 和 `comparison.svg`；不加时只写 JSON。CSV 由 C++ 计算，Python 负责检查与绘图。`results/` 保存归档快照，新结果写入 `build-bpm/`。

## 改参数观察结果

每次只改 [main.cpp](main.cpp) 中的一项，然后重新编译：

| 修改 | 预期结果 |
| --- | --- |
| `s.trailing.distance`：1.22 → 2.44 m | 有效分量和总谱降低约 6.0206 dB，来自距离平方反比 |
| `s.span`：0.4572 → 0.9144 m | 同一均匀翼段假设下升高约 3.0103 dB |
| `p.lammod`：1 → 0 | 移除层流分量，观察其贡献较大的频带 |
| 改 `s.speed` | 边界层、Strouhal 数和谱峰一起变化，不能只按 M⁵ 平移频谱 |

```sh
cmake --build build-bpm --config Release --target aeroacoustics_bpm1989 -j 4
./build-bpm/aeroacoustics_bpm1989 build-bpm/modified.csv
```

Windows 添加 `.exe`；Visual Studio 使用 `build-bpm/Release/` 路径。修改后的工况去掉 `--check`，也不要用固定工况的 `analyze.py` 校验：CSV 中的 D2 参考列仍属于原输入。比较新结果的声源列和 `total_dB`；参数扫描时为每个工况保存输入。
