# 源码目录

源码按物理模型和程序职责组织。公开接口位于仓库根目录的 `include/`，命令行程序入口位于 `apps/`。

| 目录 | 文件 | 职责 |
| --- | --- | --- |
| `acoustics/` | [empirical_models.cpp](acoustics/empirical_models.cpp) | BPM、Lowson、Guidati 等经验模型及其公式 |
| | [tno.cpp](acoustics/tno.cpp) | TNO 尾缘噪声与模型积分装配 |
| | [spectrum.cpp](acoustics/spectrum.cpp) | 参数检查、声源选择、A 计权和声能叠加 |
| | [geometry.cpp](acoustics/geometry.cpp) | 观察点与叶片前后缘坐标变换 |
| | [driver.cpp](acoustics/driver.cpp) | 声学时间状态、边界层表、湍流强度和声学输入读取 |
| `numerics/` | [quadrature.cpp](numerics/quadrature.cpp) | 61 点 Gauss–Kronrod 积分及误差估计 |
| | [backend.cpp](numerics/backend.cpp) | 普通/MKL 后端标识与向量点积 |
| `turbine/aerodynamics/` | [bem.cpp](turbine/aerodynamics/bem.cpp)、[unsteady.cpp](turbine/aerodynamics/unsteady.cpp) | BEM 诱导求解、偏斜修正和非定常翼型状态 |
| `turbine/structure/` | [blade_dynamics.cpp](turbine/structure/blade_dynamics.cpp) | 叶片模态、运动学、广义载荷与加速度 |
| `turbine/coupling/` | [mesh.cpp](turbine/coupling/mesh.cpp) | 气动和结构网格的运动、载荷传递 |
| | [rotor.cpp](turbine/coupling/rotor.cpp) | 叶素状态、气动力和转子装配 |
| | [solver.cpp](turbine/coupling/solver.cpp) | 广义 α 时间推进及结构迭代 |
| `turbine/io/` | [case_input.cpp](turbine/io/case_input.cpp) | 输入解析、配置检查、翼型表及稳态风场 |
| `interfaces/` | [c_api.cpp](interfaces/c_api.cpp) | 对外 C ABI、错误与调用状态 |
| `apps/` | [section_main.cpp](apps/section_main.cpp) | 单截面频谱示例入口 |
| | [turbine_main.cpp](apps/turbine_main.cpp) | 整机声学命令行入口和结果输出 |

`apps/turbine_main.cpp` 调用整机求解器，并将求得的节点状态交给声学驱动。整机求解器负责协调气动、结构与网格传递；声学模型调用数值积分和后端计算。

## 构建文件

- 根目录 `CMakeLists.txt`：构建开关、编译标准、MKL 查找和输出位置。
- [CMakeLists.txt](CMakeLists.txt)：声学静态库及动态库的共用源文件清单。
- [turbine/CMakeLists.txt](turbine/CMakeLists.txt)：整机模块库。
- [apps/CMakeLists.txt](apps/CMakeLists.txt)：两个可执行程序。
- 根目录 `tests/CMakeLists.txt`：数值对照探针。

重排不改变命名空间、公开头文件路径、目标名称或可执行文件位置。构建及运行命令见根目录 README。
