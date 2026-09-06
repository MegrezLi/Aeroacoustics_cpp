# 本机工具链记录

2026-09-06，按用户要求将本项目切换到 GCC / G++ / GFortran **16.2.0**。

- 发布方：[WinLibs](https://winlibs.com/)，x86_64、POSIX、SEH、UCRT。
- 发布标签：`16.2.0posix-14.0.0-ucrt-r1`，MinGW-w64 14.0.0。
- 安装目录：`D:/Code_Configuration/gcc-16.2.0/mingw64/bin`。
- 压缩包 SHA256：`c1f52294597c0b73786b2a78eb5d176d89226d2f21875eab75e783a8b1cefcc4`，已与官方 GitHub release asset 校验值核对。
- 原 GCC 7.3 保留在 `D:/Code_Configuration/mingw64/bin`。
- CMake 4.1.0，Intel oneMKL 2025.2.0。

`tools/use_gcc.ps1` 将新编译器加入当前 PowerShell 会话，并设置 CC、CXX、FC。安装没有修改系统或用户的永久 PATH。

OpenFAST `DOUBLE_PRECISION=ON` 实际启用 `-fdefault-real-8 -fdefault-double-8`。GFortran 16 已原生编译 `cotan`，源文件没有 GCC 7 的 `1/tan` 兼容改写。

首次整机尝试直接链接已有 Conda BLAS/LAPACK DLL，启动时报 Windows 32 位伪重定位越界。改为链接 MKL 正式导入库 `mkl_rt.lib` 后，原版整机正常完成 20 秒算例。因此正式原版与 C++ 接入版使用相同 MKL BLAS/LAPACK 配置。

`validation-gcc16-*.json` 为当前验证；`validation-portable.json` 与 `validation-mkl.json` 是更新工具链前的历史记录。
