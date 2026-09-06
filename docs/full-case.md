# 官方整机算例接入与复现

依赖：CMake、GCC/G++/GFortran、Intel oneMKL、Python 3.11+（脚本）、NumPy（对比）。以下命令在项目根目录的 Windows PowerShell 运行。

## 1. 下载固定版本

```powershell
python tools/fetch_openfast.py
. ./tools/use_gcc.ps1
$mklRoot = 'C:/Program Files (x86)/Intel/oneAPI/mkl/2025.2'
$env:PATH += ";$mklRoot/bin"
$env:MKL_THREADING_LAYER = 'SEQUENTIAL'
$cppRoot = $PWD.Path.Replace('\','/')
```

下载官方案例的全部 101 个输入文件并验证 Git blob 哈希；不复制历史结果作为本次输出。OpenFAST 源码放在 `_work/openfast`，算例放在 `_work/r-test/glue-codes/openfast/IEA_LB_RWT-AeroAcoustics`。`_work` 不上传 GitHub。

## 2. 编译并运行原版 Fortran

初次下载不需要恢复操作；如果已接入过 C++，先执行 `python tools/integrate_openfast.py --restore`。此操作只恢复与本项目已知补丁完全一致的两个文件，遇到其他修改会拒绝覆盖。

```powershell
cmake -S _work/openfast -B _work/build-full -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DDOUBLE_PRECISION=ON '-DCMAKE_POLICY_VERSION_MINIMUM=3.5' "-DBLAS_LIBRARIES=$mklRoot/lib/mkl_rt.lib" "-DLAPACK_LIBRARIES=$mklRoot/lib/mkl_rt.lib" -DBLAS_FOUND=TRUE -DLAPACK_FOUND=TRUE
cmake --build _work/build-full --target openfast -j 4
New-Item -ItemType Directory -Force _work/bin
Copy-Item _work/build-full/glue-codes/openfast/openfast.exe _work/bin/openfast-fortran.exe
python tools/run_full_case.py --executable _work/bin/openfast-fortran.exe --label baseline
```

## 3. 接入 C++ 并重新运行

```powershell
python tools/integrate_openfast.py --apply
cmake -S _work/openfast -B _work/build-full "-DAEROACOUSTICS_CPP_SOURCE=$cppRoot" -DAEROACOUSTICS_USE_MKL=ON -DAEROACOUSTICS_BUILD_EXAMPLES=OFF "-DMKL_DIR=$mklRoot/lib/cmake/mkl"
cmake --build _work/build-full --target openfast -j 4
Copy-Item _work/build-full/glue-codes/openfast/openfast.exe _work/bin/openfast-cpp.exe
python tools/run_full_case.py --executable _work/bin/openfast-cpp.exe --label cpp
python tests/compare_full_case.py _work/runs/baseline _work/runs/cpp
```

运行脚本为每次运行复制输入，并保存输入哈希、可执行文件哈希、返回值和日志。已有标签目录会拒绝覆盖；重新运行请指定新标签。

## 接口边界

调用链：`OpenFAST → AeroDyn → AA_CalcOutput → 七个 ISO_C_BINDING 包装例程 → aeroacoustics_kernel → C++ 公式`。

替换的例程为 `LBLVS`、`TBLTE`、`TIPNOIS`、`InflowNoise`、`BLUNT`、`Simple_Guidati`、`TBLTE_TNO`。BPM 厚度与曲线等辅助公式也由 C++ 内核内部调用。TNO 在 C++ 内部执行 MKL 向量指数和 BLAS 积分求和。

Fortran 包装例程只打包参数及复制结果。OpenFAST 外层的输入解析、状态、气动／结构求解、节点／观察点循环、几何和输出装配仍保留，确保此算例使用真实整机状态。相应辅助声学算法另有独立 C++ 实现，可通过 `AcousticDriver` 使用，但本补丁没有把 OpenFAST 全求解器替换为 C++。

日志末尾 `C++ acoustic kernel calls: ...` 记录实际调用计数。原官方案例选择 BPM、入流 Guidati 和钝度噪声；TNO、层流与叶尖开关的其他分支通过独立数值测试覆盖，不能把该默认案例的运行说成 TNO 整机验证。

## 对比标准

- 输入 101 个文件逐个哈希一致。
- 四个声学文件的通道、单位、全部 201 个时间点一致，覆盖 0–20 秒。
- 比较 OASPL、总频谱、分声源频谱及分叶片节点声级，容差 2e-5 dB（文本输出精度）。
- 二进制 `.outb` 仅归一化说明文字中的运行生成时间戳，其余字节（含全部动力学数据）要求完全一致。
- C++ 调用计数必须大于零，两个运行都必须正常结束。

测试数值体还使用独立的 Fortran 双精度参考及 Python 参考，结果见 `validation-gcc16-*.json`。这里的正确性指移植与上游公式一致，不是对风机实测噪声的独立验证。

## 2026-09-06 实测结果

当前本机运行目录为 `_work/runs/fortran-gcc16-mkl` 和 `_work/runs/cpp-gcc16-mkl`，相应可执行文件保存在 `_work/bin`。四个声学文件分别比较 402、13,668、95,676、36,180 个值，总计 145,926 个值，文本数值差异均为零。C++ 实际调用 106,128 次。

两个 `.outb` 文件均为 122,547 字节，差异仅在说明文字的生成时刻；归一化该字段后 SHA256 相同。实测结果及哈希见 [validation-full-case.json](validation-full-case.json)。

原版／接入版单次墙钟时间约 10.07／10.27 秒，运行时系统负载不同，不能据此推断性能提升或退化。本工作验证数值移植；MKL 对 TNO 的性能价值需要另做代表性基准测试。
