# CurveFitter — Qt/C++

将 ImageJ `ij.measure.CurveFitter` 的数值拟合能力移植为 C++17 / Qt Core 库，无 GUI、Java 或第三方数学库运行时依赖。

## 构建

需要 CMake ≥ 3.16、C++17 编译器和 Qt 5/6 Core 开发包。例如 Debian/Ubuntu 安装 `cmake qtbase5-dev` 或 `cmake qt6-base-dev`。

```sh
cmake -S . -B build/Qt_5_15_15-Debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/Qt_5_15_15-Debug -j2
ctest --test-dir build/Qt_5_15_15-Debug --output-on-failure
./build/Qt_5_15_15-Debug/curvefitter_example
```

在当前云环境中，开发包安装于工作区，配置命令需额外指定：

```sh
/workspace/tools/python/cmake/data/bin/cmake -S . -B build/Qt_5_15_15-Debug \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_PREFIX_PATH=/workspace/tools/qt5/usr
```

同目录下提供 `cmake` 和 `ctest`。已验证 Qt 5.15.15 / GCC 14；CMake 支持发现 Qt 6，但本次未在 Qt 6 下验证。

## 使用

```cpp
#include "curvefitter.h"

CurveFitter fit({0, 1, 2, 3}, {2, 5, 8, 11});
fit.doFit(CurveFitter::STRAIGHT_LINE);
if (fit.getStatus() == CurveFitter::SUCCESS) {
    auto p = fit.getParams(); // a=2, b=3；最后一项为目标函数 SSE
    double prediction = fit.f(p, 4); // 14
}
```

集成到其他 CMake 项目：`add_subdirectory(curvefitter)`，然后 `target_link_libraries(your_target PRIVATE curvefitter)`。

- 保留全部 25 种模型编号、公式及参数顺序：直线，2–8 次多项式，指数、幂、对数，Rodbard/逆 Rodbard，Gamma Variate，高斯，指数恢复，Chapman-Richards，误差函数。
- `setWeights()` 接收与原始输入等长的有限非负权重；零权重点不参与拟合。构造时过滤含 NaN 的点对，权重按原始索引对应；无穷值会使拟合失败。
- 非线性模型使用自动初值；复杂数据可用 `setInitialParameters()`、`setInitialParamVariations()` 设置初值与正步长。每次拟合继续使用这些设置，切换参数个数时需清空或重新设置。
- `doCustomFit()` 接受 `std::function<double(const QVector<double>&, double)>`、参数个数、说明文本、完整初值及可选步长。回调应无副作用，不能修改拟合器；异常向调用者传播。
- 拟合后读取状态，再使用参数、残差、SSE、R²、拟合优度、残差标准差及文本报告。`MAX_ITERATIONS` 时保留找到的有限最优解，调用者决定是否接受。
- 输入长度、模型编号及配置错误抛出 `std::invalid_argument`；空数据、点数不足、全零权重、奇异数据或模型定义域不兼容返回 `INITIALIZATION_FAILURE`。

## 与 ImageJ 的对应及差异

这是一份数值能力移植，不是 Java 类的逐行翻译或二进制兼容实现。

- 多项式用中心化、缩放后的 Householder QR 求解，再转换回原坐标的系数；非线性模型使用确定性的 Nelder–Mead，并解析消去最多两个线性参数。未移植 ImageJ `Minimizer` 的多线程、随机重启、ESC 中止和状态编号，收敛路径、迭代数及局部极值可能不同。每个实例同步运行，独立实例可并行使用。
- `EXP_REGRESSION` / `POWER_REGRESSION` 在对数空间拟合；允许同号的负 y。幂回归忽略 `(0,0)` 并要求正指数。`RODBARD2` 按 NIH Image 语义交换 x/y 拟合。上述三类的 SSE、R² 和拟合优度针对实际拟合空间，残差针对原始数据。
- `GAUSSIAN` 的参数 `b` 是峰值，`a` 是基线；不是直接把 `b` 当振幅。宽度归一化为正，ERF 的负宽度通过翻转幅度归一化。
- 未提供 ImageJ 宏语言字符串执行、设置弹窗、Plot 对象及 ImageJ 插件接口。自定义模型使用 C++ 回调，绘图可直接使用原始数据与 `f()` 采样值。
- `getParams()` 返回参数及末尾 SSE；本库状态枚举仅含成功、初始化失败、达到迭代上限。R² 的常量数据约定与 ImageJ 一致，返回 0。
- 本库至少要求有效点数等于参数数；不接受欠定拟合。零权重点在计算前过滤。重复 x 或不可辨识参数应由调用者检查；本库不提供置信区间或参数协方差。

## 验证与来源

CTest 运行一个测试程序，覆盖 25 种模型的默认/显式初值拟合、预测恢复、NaN 过滤、零权重、负 y 对数回归、幂回归原点、非法定义域、自定义函数、奇异直线及迭代上限。带噪指数参数另与 SciPy `least_squares` 独立参考值比较。当前共 1567 个检查通过；尚未进行 Java ImageJ 运行时的逐项差分测试，也不承诺任意数据均收敛到全局最优。

上游：[ImageJ CurveFitter.java](https://github.com/imagej/ImageJ/blob/master/ij/measure/CurveFitter.java)，本次读取文件的 SHA-256：`b71fc14eb53560bba71e8921653b0dcbafaf9f06d43edf4e7cd35f571365aced`。
保留了上游公有领域声明，见 [imagej-license.txt](imagej-license.txt)。
