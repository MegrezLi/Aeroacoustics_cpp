# 整机模块与支持范围

`aeroacoustics_turbine` 的运行入口是 `src/apps/turbine_main.cpp`。CMake 仅启用 CXX，默认构建不链接 OpenFAST、Fortran、BLAS 或 MKL。输入来自官方算例文件，风机状态从初始条件逐步计算。

表内 `.cpp` 路径相对于 `src/`。

| C++ 文件 | 对应模块与计算内容 |
| --- | --- |
| `turbine/io/case_input.cpp` | 输入字段、翼型与坐标文件、AirfoilInfo 线性插值、InflowWind 恒定风 |
| `turbine/structure/blade_dynamics.cpp` | ElastoDyn 叶片三模态方程、运动学、广义质量与载荷 |
| `turbine/aerodynamics/bem.cpp` | BEMTUncoupled 诱导公式、工作区间及入流角求根、偏斜修正 |
| `turbine/aerodynamics/unsteady.cpp` | UA_Mod=3 的附着流、分离、涡升力、滤波及离散状态 |
| `turbine/coupling/mesh.cpp`、`include/turbine/math.hpp` | NWTC 运动与载荷映射、旋转插值及坐标变换 |
| `turbine/coupling/rotor.cpp` | 转子盘与叶素坐标、风速、气动力、UA/BEM 状态传递 |
| `turbine/coupling/solver.cpp` | 气动／结构调用顺序和完整状态管理 |
| `turbine/coupling/integrator.cpp`、`structural_adapter.cpp` | 通用广义 α 积分、Newton 及叶片结构适配 |
| `turbine/coupling/modules.cpp` | 支持的模块组合、能力声明与自由度布局 |
| `turbine/simulation/` | 仿真调度、声学节点适配、声能聚合及结果输出 |
| `apps/turbine_main.cpp` | 命令行参数解析和库入口调用 |

## 官方配置

3 个叶片，叶片柔性长度 63 m，轮毂半径 2 m。每片有 17 个结构积分节点，加上根尖运动节点；气动网格为 30 个预弯节点。每片启用两个挥舞自由度和一个摆振自由度，共 9 个自由度。

结构方程包括模态形状及结构扭角、轴向缩短、质量矩阵、刚度、阻尼、重力、离心力和科里奥利项。转子以 10.04 rpm 匀速旋转，固定桨距 1.17°、预锥 −2.5°、轴倾 −5°。塔架和其余自由度锁定，与官方输入一致。

气动使用 `BEM_Mod=1`、`Wake_Mod=1`、`DBEMT_Mod=0`、`UA_Mod=3`、`FLookup=True`、`AFTabMod=1`。上游固定版本的 AirfoilInfo 默认插值阶数为 **1**，输入注释中的“默认 3”不适用于该版本。

## 每步求解顺序

1. 依据广义 α 格式预测叶片模态位移、速度及方位角。
2. 将预测运动映射至气动网格，用上一时刻输入推进非定常翼型状态。
3. 计算当前 BEM、气动力与声学输入。BEM 保留上一步的入流角，按原工作区间缩小括区，采用 Brent–Dekker 求根。
4. 保持该步气动载荷，在结构迭代中更新载荷作用点映射，求出结构加速度并校正位移、速度。
5. 按 `DT_AA` 输出声学结果；结构状态按每个 `DT` 输出。

这种调用顺序对应原 FAST_Solver 的 AeroDyn Option 2 和紧耦合 ElastoDyn。初始物理及算法加速度均设为零，与原求解器 Step0 后的状态一致。

默认参考模式下，C++ 将固定塔架条件下的结构方程分解为三个 3×3 系统，每次 Newton 迭代重新计算数值 Jacobian，修正量阈值为 1e−9，最多 12 次。原 FAST_Solver 使用包含模块输入的整体 Jacobian、缓存更新及不同的停止准则。因此 `ConvTol`、`MaxConvIter`、`DT_UJac`、`UJacSclFact` 不控制本实现的内部 Newton 迭代。`RhoInf`、`DT` 参与广义 α 系数计算；要求 `ModCoupling=3`、`NumCrctn=0`。

可通过 `SolverOptions` 修改内部迭代参数，或启用 `scaled` 模式的尺度化差分、残差判停和同一步内 Jacobian 复用；详见 [P7 求解说明](performance-p4-p7.md)。

每次 Jacobian 扰动只改变一片叶片，使用 `Rotor::structural_acceleration()` 共享该状态的节点运动，更新载荷映射并计算加速度。每轮 Newton 迭代的单叶片映射次数从 30 次降为 12 次，积分公式、扰动步长和停止准则不变。偏斜模型开关及系数在 `Rotor` 初始化时缓存。

声学驱动通过 `step_blocks()` 回调逐块交付观察点频谱，由 `AcousticWorkspace` 持有当前块，`AcousticAggregator::append()` 累计请求的输出。完整 `step_view()` 接口仍保留。工作区先准备各节点的边界层和声源谱形，再分别计算观察点的距离与指向性；主程序汇总输出后才进入下一步。频谱、TNO 积分和输出数组反复使用，湍流强度仍按原来的每步更新顺序推进。详见 [工作流](../aeroacoustics工作流.md) 和 [优化记录](optimization.md)。

BEM 使用 `IndToler` 和 `MaxIter`，默认双精度残差阈值 5e−10，并使用 1e−6 rad 的括区停止阈值。其插值和停止路径与原 Fortran Brent 例程不同，未复制原 `mod_root1dim.f90`。这些数值实现差异在整机误差报告中保留。

## 输入限制

实现覆盖上面的官方配置及已验证的 9 m/s 风速变体。支持固定转速、桨距、稳态风及表格参数的输入读取，但两项回归通过不代表任意参数组合都已验证。

主要限制由 `Case::validate_scope` 及各模块构造函数检查：

模块与自由度组合校验集中在 `configure_modules()`，Case 和冻结后的 TurbineModel 均执行检查。数值积分器已支持不同尺寸的耦合块，但当前整机物理后端仍仅实现下列范围，详见 [S4–S5](semantics-modules.md)。

- 只支持单转子、三叶片、相同结构与气动叶片文件；叶片三个模态必须全部启用。
- 仅支持 WindType=1、单张翼型表和线性翼型插值。
- 塔架、平台、传动链、变桨、偏航动力学及控制器不启用；固定平台偏移、初始叶片挠度、叶尖附加质量要求为零。
- 不支持动态入流、自由涡尾流、其他 UA 模型、塔影、塔架／机舱／尾翼气动力或扇区平均。
- 气动步长与主步长相同，声学步长须为其整数倍；整机驱动支持 `TICalcMeth=1`。
- 不提供 OpenFAST 的线性化、稳态求解、重启文件或 `.outb` 输出；C++ 库提供独立的进程内完整检查点，见 [仿真接口](simulation-api.md)。`OutList`、`DT_Out` 等原动力学输出设置不用于选择 C++ 的 `dynamics.csv` 通道。

声学库还含 TNO、层流与叶尖模型以及其他辅助算法，默认整机案例的开关并不覆盖所有声学分支。源码级 Fortran 对照单独验证这些公式。
