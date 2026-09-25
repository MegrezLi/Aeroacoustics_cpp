# 工程模型：非稳态运行与室外传播

[风场与闭环运行](#e1风场与闭环运行) · [传播](#e2传播接口) · [表面状态数据](#surface-data) · [接收时间与统计](#receiver-metrics)

新增功能由独立 C++ 模块执行：时间与空间变化的三分量风场、两质量传动链、发电机转矩响应、闭环变桨与偏航，以及空气吸收、地面反射和屏障衍射。原官方算例的固定转速路径仍保留。

## 运行示例

```sh
./build/aeroacoustics_turbine examples/engineering/closed_loop.fst build/engineering \
  --wind-grid=examples/engineering/gust_veer.wind \
  --controller=examples/engineering/controller.dat \
  --propagation=examples/engineering/propagation.dat
```

Windows 使用 `.exe`，将命令写为一行即可。三项开关可以独立使用；不指定时分别使用原稳态风、固定转速/桨距/偏航和参考自由场传播。相对路径按命令当前目录解析。批量库接口通过每个工况的 `RunOptions` 配置这些模块。

示例步长为 0.0015625 s，共运行 20 s，声学采样间隔仍为 0.1 s。`closed_loop.fst` 复用官方叶片、翼型和声学文件，源模型声速为 343 m/s，与传播示例一致。`gust_veer.wind` 是便于复现的确定性阵风与风向随高度变化数据，不是经过湍流谱和空间相干性校准的 TurbSim 场。控制器参数用于演示闭环行为，不能直接视为该机型的设计参数。

原 `.out`、`.mask`、`dynamics.csv` 继续输出；开启控制器后新增 `operation.csv`，记录转速、方位角、传动链扭转、转矩、电功率、桨距、偏航及其导数。`run.json` 中的 `engineering` 保存所用模型和控制/传播参数。`aerodynamic_torque_Nm` 是本时间步传动链所使用的上一气动时刻载荷。

## E1：风场与闭环运行

### 风场

[WindField](../include/turbine/wind.hpp) 通过 `at(time, position)` 返回全局坐标系中的 `[u,v,w]`，单位 m/s。当前 `GridWind` 读取完整的时间/x/y/z 张量网格，使用多线性插值。可接入测量、CFD 或湍流生成器导出的三分量数据，变化的横向分量能够表示风向随高度变化。求值无随机采样，同一输入可复现。

文件格式：

```text
AA_WIND_GRID 1
nt nx ny nz
t0 ... t(nt-1)
x0 ... x(nx-1)
y0 ... y(ny-1)
z0 ... z(nz-1)
u v w
... 共 nt*nx*ny*nz 行
```

时间以秒、位置以米表示。数据按时间、x、y、z 排列，z 最快；所有轴严格递增。单元素轴明确表示该方向均匀，单时间点表示定常场。其他轴越界直接报错，不端点钳制或周期回绕。风场必须覆盖全部叶片运动位置与运行时间。当前不直接读取 `.bts`，也不自动生成具有指定湍流谱/相干性的随机场；`WindType` 原文件仍保留为官方基础配置，`--wind-grid` 显式替代其查询。

风场传入 `Rotor::wind_at()`，用于各叶素相对速度、BEM、UA 和声学节点输入。Lowson 模型的湍流强度仍遵循原声学输入 `TI/AvgV` 及 `TICalcMeth=1` 约定，不根据短时间网格数据自动估计湍流统计量。用户应给出与所用风场对应的统计参数。

### 传动链与发电机

[OperatingController](../include/turbine/control.hpp) 持有两质量传动链和执行器的全部历史状态。设传动比 `N=ωg/ωr`，转子侧轴扭转 `δ=θr−θg/N`：

```text
Ts = K δ + D (ωr − ωg/N)
Jr dωr/dt = Taero − Ts
Jg dωg/dt = Ts/N − Tg
dδ/dt = ωr − ωg/N
dθr/dt = ωr
Pe = η Tg ωg
```

`Jr` 是包括叶片和轮毂在内的转子侧等效惯量，`Jg` 为高速轴侧发电机惯量；`K/D` 均折算到低速轴。齿轮箱采用理想刚性传动比，电效率作用于输出功率。发电机电磁转矩以一阶时间常数响应指令，并限制最大值和变化率。

气动转矩由各叶素力及力矩沿参考弧长积分得到。线性力与线性位置的乘积采用精确积分，和已有 `LoadMap` 的合力矩一致。积分得到的方位角用于叶片基底，变转速的欧拉惯性项、变桨和偏航产生的角速度/角加速度，以及偏航导致的轮毂平移速度/加速度，均传入结构运动学。

当前采用常数等效转子惯量，未把叶片弹性变形引起的瞬时惯量变化及其全部惯性反作用反馈至传动链；塔架/平台仍固定。发电机采用转矩响应模型，不含电磁电路、电网故障或电能质量模型。

### 控制律与执行器

- **转矩控制**：滤波后的转速用于 `Kopt ωg²` 指令，以当前额定电功率对应的转矩为上限，并受独立转矩限值约束。
- **集体变桨**：PI 根据滤波后的转子等效转速与目标转速之差调节桨距，包含积分抗饱和、角度限值和变化率限制；三片叶片使用相同指令。
- **变桨执行器**：二阶伺服响应，参数为固有角频率、阻尼比及最大速率。
- **偏航控制**：跟踪轮毂处水平风向，使用最短角差、死区和二阶速率受限伺服。它是执行器响应模型，不是塔顶偏航载荷方程。
- **降噪模式**：到 `NoiseStart` 后切换额定转速比例、功率比例和最小桨距。声学结果仍由实际运动、气动和声源模型计算，未人为扣除固定 dB。降噪目标不保证所有观察点、频带和风况都降低。

各控制参数的单位见 [controller.dat](../examples/engineering/controller.dat)。这是本项目的简化闭环控制器，并非 ROSCO 或 ServoDyn 的完整移植；原 `CompServo`、`GenDOF` 等 OpenFAST 开关不能用于选择它。调用方须显式使用 `--controller` 或 `SolverOptions::controller`。

每个结构步内，控制/传动链采用 RK4 子步积分，子步同时受用户 `MaxStep`、轴系固有频率、阻尼衰减率和执行器时间尺度限制。该步气动转矩与风向保持为上一个气动时刻的值。随后重新计算当前气动，并校正叶片模态。这个分步耦合仍需要主步长收敛检查，减小内部控制步长不能消除气动/结构交换误差。

`Solver`/`Simulation` 的复制、检查点和重置均涵盖控制器积分、滤波、传动链、执行器、UA/BEM 和声学 TI 历史。求解失败后不得继续推进失败状态。转子或发电机停止/反转会报错，当前不包含启动、停机和机械制动流程。

## E2：传播接口

[OutdoorPropagation](../include/propagation.hpp) 在各节点、各观察点、各频带聚合之前处理传播。源模型已经包含 `1/r²` 和指向性，因此传播层不再给直达声重复施加几何衰减。

### 空气吸收

均匀温度、相对湿度和气压下，按 ISO 9613-1 的氧/氮弛豫解析关系计算中心频率吸收系数 `α(f)`，单位 dB/m，然后对实际源边缘到受声点的路径扣除 `αr`。环境参数要求 −20 至 50 °C、10% 至 100% RH 和 50–120 kPa。

原算例 34 个频带覆盖 10–20000 Hz；其中低于 50 Hz、高于 10000 Hz 的频带属于该关系在标准主要范围之外的延伸计算。这里使用中心频率近似，尚未按带内源谱对吸收系数重新积分。传播声速控制相位，源模型声速控制马赫数；均匀环境中应使用一致的值。

### 地面

支持无反射、刚性平面和给定归一化复阻抗三种条件。地面标高为全局 `GroundZ`，源和受声点不能低于该平面。反射路径通过镜像受声点重新调用原声源模型，保留其实际路径长度和发射方向；没有复用直达声的指向性。

阻抗输入为 `Z/(ρc)`，采用 `exp(+iωt)` 约定，要求实部为正。平面波反射系数为 `(Z cosθ−1)/(Z cosθ+1)`。当前复阻抗是频率无关的理想边界；并未实现真实土壤的频率相关阻抗或球面波地面修正。

直达与反射声按同一声源的相干路径叠加；对频带内相位差解析积分，交叉项含 `sinc(π Δf Δr/c)`，避免仅在中心频率计算干涉产生频带混叠。它假设带内源密度及反射幅值近似恒定，并采用相同源相位；经验声源模型本身不提供完整复数辐射场。不同节点和机制继续非相干能量叠加。

### 地形遮挡

`NumScreens` 后可给出山脊/薄屏障的水平顶边端点及顶标高，格式 `x1 y1 x2 y2 top_z`。传播路径与其水平投影相交时，根据相对直达线的高度和声波波长计算单刃 Fresnel 参数，使用标量刀刃衍射近似；多个屏障取主导衰减，不将各项直接相加。直达路径和经过镜面反射点的两段反射路径分别检查遮挡，反射路径同样取主导衰减。地面/屏障联合计算仅修正各路径幅值，相位仍按镜像路径近似，未加入绕射相位。

这种表示适用于以孤立主导山脊或薄屏障近似的路径，不等于任意数字地形模型，也不计算有限屏障侧绕射、多重绕射或建筑群反射。未提供气象梯度折射、声影区射线或抛物方程求解；本版本按均匀介质直线路径工作，不能宣称完整实现 ISO 9613-2:2024。

## 验证与使用边界

- 解析测试：四维线性风场、两质量传动链的非对角耦合响应、气动转矩与载荷映射守恒、运动位置的速度/加速度差分。
- 传播测试：无损刚性地面在地面受声点增加 6.0206 dB、路径吸收、遮挡衰减、非法阻抗/频带和输入范围。
- 整机测试：正常/降噪模式、控制状态约束、检查点重放、四档主步长、非稳态风场越界拒绝。
- 原 Fortran/C++ 声学和整机对照继续运行；新增工程模型尚无同配置 Fortran 或现场测量对照。它们的解析与收敛测试不能替代实测验证。

结果与实际步长误差保存在 [validation-engineering.json](validation-engineering.json)。工程算例不能仅根据是否成功运行判断时间步足够小；控制器增益、惯量、传动链刚度、气动状态及风场时间分辨率都影响步长要求。

## 模型资料

两质量传动链的数据组织参考 [NREL：Dynamic Models for Wind Turbines and Wind Power Plants](https://nrel.gov/docs/fy12osti/52780.pdf)。转矩、变桨、执行器和运行降额的工程背景见 [ROSCO 参数说明](https://rosco.readthedocs.io/en/latest/source/rosco_toolbox.html) 与 [功率控制示例](https://rosco.readthedocs.io/en/latest/source/examples.html)。风场坐标和输入组织参见 [OpenFAST InflowWind](https://openfast.readthedocs.io/en/dev/source/user/inflowwind/driver.html)。

空气吸收范围见 [ISO 9613-1](https://www.iso.org/standard/17426.html)；计算关系另可核对 [Acoustic Toolbox 的公式文档](https://acoustic-toolbox.readthedocs.io/en/latest/standards/iso_9613_1_1993/)。刀刃标量波近似见 [ITU-R P.526-14，单刃衍射](https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.526-14-201801-I!!PDF-E.pdf)，此处使用声速计算声波波长，不采用无线电地球曲率模型。ISO 9613-2 工程方法的完整适用范围见 [ISO 官方说明](https://committee.iso.org/standard/74047.html)。

<a id="surface-data"></a>

## E3：翼型、转捩与表面状态数据

```sh
./build/aeroacoustics_turbine examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/surface --surfaces=examples/engineering/surfaces.dat
```

`SurfaceSet` 在仿真初始化时读入并冻结数据。清单用 `NumStates` 指定数量，`StateFiles` 起逐行列出状态文件；文件中的相对路径按其自身所在目录解析。一个状态绑定一个从 1 开始的 `AirfoilID`，作用于使用该翼型的所有叶素和三片叶片。重复编号或超出当前机型的编号报错。当前不表示同一翼型在三片叶片上的不同劣化分布。

每个状态文件的字段如下，完整例子见 [surface-test.dat](../examples/engineering/surface-test.dat)。

| 字段 | 含义 |
| --- | --- |
| `State`、`Provenance`、`UncertaintyNote` | 状态名称、数据出处和不确定度的定义，均必填 |
| `BoundaryLayer` | 使用现有 OpenFAST 格式的边界层表 |
| `Polar` | 完整翼型文件，或 `none`；完整文件包括极曲线、UA 参数和坐标引用 |
| `AlphaMinDeg/AlphaMaxDeg`、`ReMin/ReMax` | 数据已验证的有效区间，必须落在所提供表格内 |
| `TransitionSuction/TransitionPressure` | 吸力面、压力面的转捩位置 `x/c`，范围 0–1 |
| `RoughnessM`、`ErosionM` | 表面状态的粗糙度/侵蚀几何标签，单位 m |
| `RelativeUncertainty` | 0–1 的相对输入不确定度声明；其统计含义由说明字段给出 |
| `TEThicknessM`、`TEAngleDeg` | 钝度模型使用的尾缘厚度和夹角，单位 m、度 |

边界层表每行依次为攻角、两侧 `Ue/Uinf`、两侧 `δ*/c`、两侧 `δ99/c`、两侧 `Cf`；每一对均为**吸力面在前、压力面在后**。Re 按旧文件格式以百万为单位，状态文件中的 Re 则是无量纲实际数值。插值后厚度乘当地弦长，转换为米。要求 `Ue/Uinf>0`、`0<δ*≤δ99`、`Cf≥0`，表格坐标递增且数据有限。

指定状态的截面使用其边界层表，其余截面保持原 `BLMod/TripMod` 配置。提供 `Polar` 时，BEM 与 UA 一并使用新翼型；仅提供边界层时，气动力学保持原翼型。有效范围检查针对实际声学查询和提供新极曲线时的实际气动状态，越界即失败，不使用默认端点钳制掩盖超范围数据。BEM 求根试探仍可查询所提供的完整极曲线。

转捩、粗糙度与侵蚀字段用于说明数据对应的表面状态，**不会自行生成经验声级修正**。它们的影响需要体现在配套的极曲线、UA 参数和边界层中。当前也不将输入不确定度自动转换成噪声置信区间；这属于 E8。TNO 截面主路径仍使用参考外缘速度比约定，E9 尚未改变。

每次运行输出 `surface_datasets.csv`，保存数据出处、有效范围、表面标签、不确定度及是否替换极曲线。随仓库的 `surface-test-bl.dat` 是插值和输入敏感性测试数据，不能作为真实侵蚀叶片的数据使用。输入格式背景见 [OpenFAST 边界层说明](https://openfast.readthedocs.io/en/main/source/user/aerodyn-aeroacoustics/App-usage.html)。

<a id="receiver-metrics"></a>

## E5、E7：接收时间、音调输入与工程统计

```sh
./build/aeroacoustics_turbine examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/metrics --surfaces=examples/engineering/surfaces.dat --metrics=examples/engineering/metrics.dat
```

两项开关可独立使用。`--metrics` 的 `ObserverFile` 可替换整机观察点列表，示例提供四个受声点。坐标与原模型一样采用全局 m；文件第一行为点数，第二行为说明，之后每行 `x y z`。这些点都由声源模型重新计算，不对已有两点结果进行空间插值。

### 配置与输出

| 配置 | 作用 |
| --- | --- |
| `Start/End` | 请求的接收时间区间，单位 s；用于剔除初始过程、选择分析时段 |
| `ReceiverDT` | 接收时间网格间隔；减小它不能补回原 `DT_AA` 丢失的变化 |
| `RetardedTime` | 是否按各节点的实际源位置建立接收时间；关闭时使用原源时间 |
| `AMWindow`、`AMMinHz/AMMaxHz` | AM 完整窗口长度及基频搜索范围；窗口须至少覆盖最低频率的两周期 |
| `WindBinWidth` | 水平轮毂风速分箱宽度，m/s；区间为 `[kΔU,(k+1)ΔU)` |
| `BackgroundLAeq` | `none` 或给定的恒定 A 计权背景声级，与预测风机声能相加 |
| `ApparentPower` | 是否给出自由场条件下的方向性等效表观声功率 |
| `NumTones`、`ToneFiles` | 独立窄带线声源数量及输入文件；0 表示不添加 |
| `MaxValues` | 历史数据数值个数上限；超限报错，不静默截断 |

| 输出 | 内容 |
| --- | --- |
| `receiver_history.csv` | 接收时间、轮毂风速、风机 A 声级、含背景总声级和频带结果 |
| `receiver_map.csv` | 每个坐标点的实际统计时段、LAeq、L5/L50/L95 和可选表观声功率 |
| `wind_bins.csv` | 每个观察点、风速区间的有效持续时间和 LAeq |
| `am_windows.csv` | 完整窗口的主调制频率、重构峰峰值及重构 P95−P5 |
| `receiver_tones.csv` | 各独立音调的接收频率、未计权声级和相对所在宽带频带的能量比 |
| `metrics.json` | 时间基准、频带计权、背景、统计口径、简化假设及音调出处 |

原 `.out` 文件继续表示源时刻的七类叶片声源；独立输入音调加入新增接收结果，未塞入原七个通道。`NrOutFile=1` 时也能输出完整统计。CSV 用 `-inf` 表示零声能；音调所在频带没有宽带能量时，能量比留空。

### 接收时间与运动音调

`ArrivalSeries` 使用 `t_receive=t_emit+r(t_emit)/c`，要求接收时间严格递增。宽带各节点的前缘入流贡献与尾缘贡献分开延迟，随后在共同接收时间网格上线性插值**均方声压**并累加。统计只覆盖所有声源、观察点和轮毂风速历史的公共有效区间；不把缺失历史补零，也不向历史之外外推。请求区间与实际区间分别写入元数据。

宽带源模型原有的马赫数与指向性修正保留，本层不再叠加一套宽带多普勒频率或幅值修正。这是运动声源**声级包络**的延迟处理，不是复声压波形、完整可压缩运动声源求解器或音频合成。启用时延时暂不接受相干地面反射，因为两条路径需要各自的历史和相位；直接路径可以同时使用均匀空气吸收和原屏障幅值近似。源模型与传播层声速须一致。

音调文件指定 `Name/Provenance`、从 1 开始的 `Blade/Node`、`FrequencyHz`、`RotorOrder`、`ReferenceSPL` 和 `ReferenceDistance`。它附着于气动节点中心，发射频率为 `FrequencyHz+RotorOrder×ω/(2π)`。输入声级是参考距离处的未计权均方声压级，采用各向同性 `1/r²` 能量衰减；它是外部给定的独立线声源，不由 BPM 频带结果反推。

接收频率乘以相邻发射/接收时间增量比 `Δt_emit/Δt_receive`，因此恒速接近和静止极限均可直接核对。幅值按给定离散音调的声压包络输运，不另加运动单极子的幅值模型。音调按接收频率归入当前频带并单独施加 A 计权，彼此按不相干能量叠加；目前无相位、相干拍频或齿轮/电磁音调生成模型。音调的空气吸收使用接收频率和插值后的源—接收点距离计算；音调路径仅支持自由场和空气吸收，不支持地面/屏障组合。源线须与宽带模型代表的机制区分，避免人为重复输入同一声源。

### 平均、分箱、AM 与表观声功率

LAeq 对 A 计权均方声压作梯形时间积分，再除以有效时长并取对数。已有 A 计权输入不重复计权；未计权频带按原中心频率公式处理，音调使用其实际接收频率。L5/L50/L95 是按时间持续量加权的超越声级，静音不替换成 0 dB。背景是另外加入预测结果的恒定声能，未实现实测总声级的背景扣除或有效性评级。

轮毂风速按初始轮毂位置到观察点的时延对齐，作为该观察点的工况标签；不是所有分布声源的唯一发射时刻。每个时间间隔在跨越风速箱边界处切分，并积分线性声能，避免按样本数量平均。该风速未折算到 10 m 高度，也未作 IEC 标准化风速处理。

AM 对 A 声级序列作离散傅里叶分析，在设定范围内选取最强基频，重构基频及不超过 Nyquist 的二、三次谐波，输出峰峰值和 P95−P5。仅计算完整窗口；静音/不足样本的窗口标记为未解析。这些是描述性调制指标，不包含 IOA 方法的全部频段选择、趋势处理、显著性筛选和长期评级，不能称为 IOA/IEC 合规 AM。音调的“线声能/所在宽带频带声能”也不是临界带音调可听度或罚值。IOA 的正式方法见 [原报告](https://www.ioa.org.uk/sites/default/files/AMWG%20Final%20Report-09-08-2016_1.pdf)。

`ApparentPower=true` 仅允许原自由场预测，计算 `LAeq+10log10(4πr²)`，r 为初始轮毂到受声点距离。它是指定方向上的等效表观值，不是将各方向积分得到的总辐射声功率，也不是地面测量板条件下的标准值。启用传播模型时拒绝此换算，避免把已受地面/空气/屏障影响的声压直接命名为声功率。库函数另可显式输入已知反射修正，标准测量流程仍需独立实现。

[IEC 61400-11](https://webstore.iec.ch/en/publication/5428) 的发射表征与 [IEC TS 61400-11-2:2024](https://webstore.iec.ch/en/publication/62414) 的受声点测量用途不同；本模块记录计算条件，没有宣称通过这些测量标准的完整符合性验证。

### 验证与后续范围

解析检查覆盖非均匀时间步的声能平均、持续时间百分位、已知正弦调制、匀速运动音调多普勒、静止极限、跨风速箱积分和表观声功率几何换算。整机检查覆盖表面数据敏感性、输入范围拒绝、A 计权一致性、接收点地图、音调频率变化、闭环风场/空气吸收组合、检查点/重置及并发。结果见 [validation-metrics.json](validation-metrics.json)，命令见 [验证文档](validation.md#surface-metrics)。

由于尚无用户的真实翼型或音调数据，当前算例用于软件和输入敏感性验证。数据标定、真实侵蚀预测、完整宽带运动声源频率映射、音频合成、IEC 音调/背景处理及测量标准符合性仍需后续工作。


<a id="farm"></a>

## E6：准稳态多机尾流与风场声学

```sh
./build/aeroacoustics_farm examples/farm/farm.dat build/farm-results
./build/aeroacoustics_turbine examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst build/tower-results 20 --tower=examples/farm/tower-test.dat
```

风场输出目录必须不存在，父目录须已存在。Windows 添加 `.exe`。示例运行两台相距 800 m 的风机，并加入固定塔架与两个给定声功率源；布局、塔架尺寸和外部谱均为软件测试数据，不代表真实风场测量或机型标定。

### 配置与坐标

[farm.dat](../examples/farm/farm.dat) 的文件路径相对于所在配置文件解析，沿用“值在前、字段名在后”的格式。

| 输入 | 含义 |
| --- | --- |
| `Provenance` | 工况及数据来源，必须填写 |
| `Duration`、`StatisticsStart` | 运行终点和统计起点，单位 s；终点须为结构和声学时间步的整数倍，统计起点不早于声学起点 |
| `Wakes`、`WakeExpansion`、`MaxCt` | 尾流开关、无量纲线性扩张率 k>0、尾流 CT 上限（0<MaxCt<1） |
| `MaxValues` | 保留数值历史的标量数上限；包含上游历史、汇总频带/时序和外部源缓存，不含模型、容器和临时工作区的全部内存 |
| `ObserverFile` | 首行观察点数量、第二行说明，后续各行全局 x/y/z，单位 m；见 [receivers.dat](../examples/farm/receivers.dat) |
| `Propagation` | `none` 或已有传播配置；本风场入口仅接受自由场/空气吸收，拒绝地面和屏障 |
| `NumTurbines`、`TurbineFiles` | 1–1000 台；首个路径与标签同行，其余路径逐行列出 |
| `NumSources`、`SourceFiles` | 0–10000 个固定声源；数量非零时按同样格式提供路径 |

各机组文件包含唯一 `Name`、`.fst` 的 `Case`、全局 `X/Y` 原点，以及 `Controller/Tower/Surfaces` 文件或 `none`。XY 原点是单机坐标系的平移，不是轮毂位置；轮毂高度和偏置仍来自各自机型。场址为平面，原点 z=0。公共观察点和固定声源使用全局坐标；塔架文件使用本机坐标。

所有机组须有相同的稳态环境风速、参考高度、切变和水平风向，相同空气密度/声速、结构/声学步长、声学起点及频带。`NrOutFile>=2`。风向约定为 `(cos(PropagationDir), -sin(PropagationDir), 0)`；机组轴线水平投影须顺风对齐，轴倾角小于 10°，轮毂位置不移动。可以接入独立变速变桨控制，但不接受改变对齐方向的偏航运动；当前风场入口不读取网格风或时变风向。

### 尾流与推力反馈

程序按轮毂的沿风位置排序，每台机组独立运行 BEM、UA、结构和可选控制状态，完成后把时间、实际轮毂入流和推力系数历史冻结。下游查询同一时刻的上游历史，采用 Jensen 顶帽形尾流：尾流半径 `Rw=R+k*x`，尾流内速度亏损 `ΔU=Uhub*(1-sqrt(1-CT))*(R/Rw)^2`，上游或尾流外为零。形状关系可参照 [FLORIS 的 Jensen 实现](https://github.com/NatLabRockies/floris/blob/main/floris/core/wake_velocity/jensen.py)。本程序的叠加约定是各上游**有量纲**亏损的平方和开方，再从该点环境风扣除；不是将不同上游的无量纲亏损直接相加。

推力由当前叶片气动力沿叶展积分并投影到来流方向，`CTraw=T/(0.5*rho*pi*R²*Uhub²)`。轮毂风速包含上游尾流，不含本机 BEM 诱导或塔影。尾流使用 `clamp(CTraw,0,MaxCt)`，气动力与结构求解结果不被改写。启动瞬态或高推力可能超出该尾流公式的适用域，必须查看 `wake_diagnostics.csv` 的原值、使用值和限幅标志；`farm.json` 汇总限幅次数。示例默认 `MaxCt=0.95` 是显式数值设定，不是所有机型通用的标定值。

风机排除球不能相交；尾流可能覆盖下游转子的机组对，其沿风距离至少为上游直径的两倍。这只是拒绝明显近场配置的下限，不能保证任意布局精度。亏损耗尽局部入流时直接报错，不设置隐含风速下限。各叶素点分别判断是否位于尾流内，未实现平滑尾流边界或专门的重叠面积积分。

此版本不含尾流附加湍流、动态摆动、输运时延、偏航偏转、阻塞/上游反馈和地形；也没有把它们折算为经验噪声增量。它不等同于 [FAST.Farm 动态尾流模型](https://openfast.readthedocs.io/en/main/source/user/fast.farm/FFarmTheory.html)。

### 固定塔架入流影响

[tower-test.dat](../examples/farm/tower-test.dat) 包含 `Provenance`、本机坐标 `X/Y`、`Potential/Shadow` 开关、`NumStations` 和 `0 Stations` 后的 `z diameter Cd` 三列表。至少两个按 z 递增的站位，长度单位 m，直径>0，Cd 在 [0,3]。直径和 Cd 沿塔高线性插值。

势流采用圆柱公式，以塔半径归一化风向/横向距离 x、y：`du=(y²-x²)/(x²+y²)²`、`dv=-2xy/(x²+y²)²`。Powles 塔影只在下风向施加余弦平方亏损，其最大幅度限制为局部水平风速的 0.5。有限塔端在一倍半径内衰减；离塔面超过 40 倍半径或小于等于 0.02 倍半径时不施加修正，穿入排除包络时报错。公式和截断依据 [AeroDyn 塔架影响理论](https://openfast.readthedocs.io/en/dev/source/user/aerodyn/theory.html)，未实现其全部选项。

这里用**叶素位置查询到的未受本塔影响的水平风**定义局部幅值与方向，保持竖向风分量；没有另外在塔节点查询入流。修正进入叶素 BEM/UA 及后续声源计算。塔架本身没有弹性自由度、阻力载荷、辐射声或声散射，亦无机舱气动。原 `.fst/AeroDyn` 塔架选项仍受原有检查约束，新增功能通过 `--tower` 或风场机组 `Tower` 文件接入。

### 机械与冷却声源输入

[cooling-test.dat](../examples/farm/cooling-test.dat) 与 [tone-test.dat](../examples/farm/tone-test.dat) 提供固定声源接口。文件含 `Name/Provenance`、`Kind`（`bands` 或 `tones`）、全局 `X/Y/Z`、单位方向向量 `AxisX/AxisY/AxisZ`、`Directivity`、`NumFrequencies`，以及 `0 Spectrum` 后的频率/声功率两列表。频率须为正且递增。

声级是**未计权声功率 Lw，参考 1 pW**；bands 为对应频带总功率，tones 为离散谱线总功率，不是 SPL 或 PSD。bands 的频率须等于计算频带中心；tones 可在覆盖范围内任意取值，按频带边界归属。只叠加输入表中指定的贡献，缺失频带按零声能处理，不能把不完整测量谱当成完整机械总声功率。

指向性 `Q=1+beta*cos(theta)`，beta=`Directivity`∈[-1,1]，球面积分为 4π。自由场远场采用 `p²=rho*c*W*Q/(4πr²)`，声压参考 20 μPa；再按路径加入可选空气吸收。频带用中心频率、音调用精确谱线频率计算吸收和 A 计权。接收点与声源不能重合，使用者还需保证测点位于点声源远场适用范围。

这些源在一次运行内位置、频率和声功率均不变；不从齿轮、发电机、电磁激励或结构振动预测声级，也未按转速/负载插值。应输入对应工况的实测或可信预测谱，并避免与已有叶片机制重复。当前示例均为合成数据。

### 汇总与输出

每台机组的普通输出保存在 `turbine_N/`，N 对应输入清单序号；求解顺序由来流位置决定。汇总使用共同**源时间**的 A 计权频带均方声压，已经 A 计权的输入不重复计权。所有机组和固定源按不相干声能相加，没有相位或交叉项。

| 输出 | 内容 |
| --- | --- |
| `farm_history.csv` | 源时间、观察点、总 A 声级及各 A 计权频带 |
| `farm_receivers.csv` | 全局坐标、实际源时间窗、LAeq 和时间加权 L5/L50/L95 |
| `source_contributions.csv` | 每台风机/固定源在各观察点的 LAeq，便于按声能核对总量 |
| `wake_diagnostics.csv` | 每个结构步的推力、轮毂风速、原始/使用 CT 与限幅标志 |
| `stationary_sources.csv` | 固定源位置、轴向、指向性、输入频谱及来源 |
| `farm.json` | 尾流设定、限幅次数、介质、时间窗、时基与叠加假设、机组顺序；全部输出成功后写入 |

统计对窗口端点的声能插值，再作时间积分；不会平均 dB。这里未调用 E5 接收时延重采样，因此风场时序不能用于完整的接收时间 AM、相干音调或音频评价。固定源和各机组的 LAeq 能量可相加，L5/L50/L95 不可直接相加。真实工程使用还需要输入谱、尾流扩张率与推力适用域的标定，以及足够长的统计窗口；当前测试不证明实测精度。验证证据见 [E6 验证](validation.md#farm)。
