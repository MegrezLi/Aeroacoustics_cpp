# 工具链

默认独立 C++ 程序需要 CMake ≥3.20 和 C++17 编译器。GFortran 只用于原源码对照；MKL 为可选声学后端。

## 本机 GCC

已安装 GCC / G++ / GFortran 16.2.0，来自 [WinLibs](https://winlibs.com/)，x86_64、POSIX、SEH、UCRT，MinGW-w64 14.0.0。

- 发布标签：`16.2.0posix-14.0.0-ucrt-r1`。
- 安装位置：`D:/Code_Configuration/gcc-16.2.0/mingw64/bin`。
- 安装包 SHA256：`c1f52294597c0b73786b2a78eb5d176d89226d2f21875eab75e783a8b1cefcc4`，已与发布文件校验值核对。
- 原 GCC 7.3 保留在 `D:/Code_Configuration/mingw64/bin`。
- 验证使用 CMake 4.1.0；原 Fortran 的 BLAS/LAPACK 使用 Intel oneMKL 2025.2。

`. ./tools/use_gcc.ps1` 设置当前终端的编译器及 PATH，不修改永久系统配置。只有执行 Fortran 对照时才需要 GFortran。

## 可选 MKL 声学后端

```powershell
. ./tools/use_gcc.ps1
$mklRoot = 'C:/Program Files (x86)/Intel/oneAPI/mkl/2025.2'
$env:PATH += ";$mklRoot/bin"
$env:AEROACOUSTICS_RUNTIME_DIRS += ";$mklRoot/bin"
$env:MKL_THREADING_LAYER = 'SEQUENTIAL'
cmake -S . -B build-mkl -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DAEROACOUSTICS_USE_MKL=ON "-DMKL_DIR=$mklRoot/lib/cmake/mkl"
cmake --build build-mkl -j 4
./build-mkl/aeroacoustics_example.exe build-mkl/spectrum.csv
python tests/compare_fortran.py build-mkl/libaeroacoustics_shared.dll --fortran-compiler $env:FC --report build-mkl/kernels.json
```

MinGW 通过 `mkl_rt` 单一动态接口链接；运行时应能找到 MKL DLL。TNO 积分调用 `vdExp`、`cblas_dgemv` 和 `cblas_ddot`，是否加速取决于工况和规模。

## 原 Fortran 整机参考

下载固定原源码后，使用同一 PowerShell 会话：

```powershell
cmake -S _work/openfast -B _work/build-reference -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DDOUBLE_PRECISION=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 "-DBLAS_LIBRARIES=$mklRoot/lib/mkl_rt.lib" "-DLAPACK_LIBRARIES=$mklRoot/lib/mkl_rt.lib" -DBLAS_FOUND=TRUE -DLAPACK_FOUND=TRUE
cmake --build _work/build-reference --target openfast -j 4
python tools/run_full_case.py --executable _work/build-reference/glue-codes/openfast/openfast.exe --label fortran-reference
```

OpenFAST 双精度配置使用 `-fdefault-real-8 -fdefault-double-8`。原 Fortran 参考链接 MKL；独立 C++ 默认构建不链接 MKL，也不复用这个 OpenFAST 可执行文件。两者的差异保存在验证报告中。
