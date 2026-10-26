# 运行快照

这组结果于 2026-10-08 在 Windows、GCC 16.2.0、CMake Release 配置下生成。
归档声学库基于提交 `90bc40778fed432061d449177550e5273b35bfb9`。
CSV 由同一算例中的 `main.cpp` 计算，JSON 与图像由 `analyze.py` 生成。

| 文件 | 用途 |
| --- | --- |
| `comparison.csv` | 27 个频带的各声源、总谱、Table D2 总声级和差值 |
| `summary.json` | 从 CSV 核算的误差与检查结果 |
| `comparison.png` | 教程中使用的对照图 |
| `comparison.svg` | 同一图的矢量版本 |

Table D2 数值来源为 NASA-RP-1218 印刷页 114；27 个总声级已与原表图像逐行核对。
这里保存的是模型输出对照。参数和重新生成命令见[算例教程](../README.md)。
