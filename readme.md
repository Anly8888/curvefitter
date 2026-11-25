# CurveFitter — Qt/C++

按 ImageJ 原版数值流程移植 `ij.measure.CurveFitter` 与 `ij.measure.Minimizer`。此前简化的 Nelder–Mead 和多项式 QR 路径已删除；全部模型统一使用原版的回归消元、单纯形及参数转换流程。

## 构建与测试

默认只需要 C++17、CMake ≥ 3.16 和 Qt 5/6 Core；单元测试使用 Qt Test。不依赖 Qt Widgets。

使用 Qt 5 时，Debian/Ubuntu 的基础开发依赖为 `cmake qtbase5-dev`。
Java 仅用于重新生成对照数据，库及常规测试不需要 Java 或联网。

```sh
cmake -S . -B build/Qt_5_15_15-Debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/Qt_5_15_15-Debug -j2
ctest --test-dir build/Qt_5_15_15-Debug --output-on-failure
```

当前云环境的开发包安装在工作区，配置时需指定：

```sh
/workspace/tools/python/cmake/data/bin/cmake -S . -B build/Qt_5_15_15-Debug \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_PREFIX_PATH=/workspace/tools/qt5/usr
```

同目录提供 `cmake` 与 `ctest`。已实测 Qt 5.15.15 / GCC 14；Qt 6 尚未实测。代码按仓库 `.clang-format` 格式化，并按 AGENTS 手动检查函数间距和逻辑分段。

### 可选：示例程序

默认不构建示例，不影响库和单元测试。通过 `-DCURVEFITTER_BUILD_EXAMPLES=ON` 启用：

```sh
cmake -S . -B build/Qt_5_15_15-Debug -DCMAKE_BUILD_TYPE=Debug \
  -DCURVEFITTER_BUILD_EXAMPLES=ON
cmake --build build/Qt_5_15_15-Debug --target curvefitter_example
./build/Qt_5_15_15-Debug/curvefitter_example
```

### 可选：字符串公式

仅在需要执行字符串公式时启用 `-DCURVEFITTER_ENABLE_STRING_FORMULAS=ON`。
启用时需要与 Qt Core 同版本的 Qt Qml 开发包（Debian/Ubuntu 上 Qt 5 对应 `qtdeclarative5-dev`）。

该选项默认关闭。关闭时不查找或链接 Qt Qml，且不提供接受字符串公式的
`doCustomFit()` 和 `getNumParams()` 重载。使用 CMake 链接 `curvefitter` 的调用方
可通过 `#ifdef CURVEFITTER_ENABLE_STRING_FORMULAS` 检查该能力。
内置模型、`Minimizer` 和 C++ 回调形式的自定义拟合不受影响。

## 拟合接口

```cpp
#include "curvefitter.h"

CurveFitter fitter({0, 1, 2, 3}, {2, 5, 8, 11});
fitter.doFit(CurveFitter::STRAIGHT_LINE);
if (fitter.getStatus() == Minimizer::SUCCESS) {
    auto p = fitter.getParams(); // a=2、b=3；末尾为原版定义的 SSE
    double predicted = fitter.f(4); // 14
}
```

CMake 集成：`add_subdirectory(curvefitter)` 后链接 `curvefitter`。支持原版全部 25 种模型编号、公式与参数顺序，提供原始数据、残差、SSE、标准差、R²、拟合优度、宏公式文本和结果报告。

- `setInitialParameters()`、`setInitialParamVariations()` 设置初值及变化尺度；省略时按原版估计。空 Qt 数组对应 Java 的 `null`。
- `setWeights()` 沿用原版权重数组和计算顺序，不剔除零权重点。权重与拟合使用的数据顺序一致；输入中的 NaN 点对在构造时过滤，建议调用者在传入权重前统一清洗数据。权重通常取 `1/error²`。
- 回调形式的 `doCustomFit()` 接受 `std::function<double(const QVector<double>&, double)>`。通过 `setOffsetMultiplySlopeParams()` 指定可由回归直接求出的偏移、乘法因子或斜率。
- `getMinimizer()` 返回实际使用的优化器引用，可设置随机种子、最大线程数、参数分辨率、重启次数及迭代上限，也可调用 `abort()`。
- 为遵循上游，每个实例保留拟合器和优化器的原始状态生命周期；独立实验或需要清零迭代计数时，创建新实例。配置与结果访问不要与正在运行的拟合并发，`abort()` 可跨线程调用。

## Minimizer

`minimizer.h` 可独立使用。目标函数接收可修改参数数组，前 n 项是搜索参数，第 n 项留给目标值；额外元素通过 `setExtraArrayElements()` 分配，供目标函数保存辅助结果。

保留了原版的：

- 与 `java.util.Random` 相同的 48 位随机数算法和 `nextFloat` / `nextDouble` 顺序；随机正交单纯形与 NaN 定义域搜索。
- 反射、扩张、内外收缩、整体缩小、中心点检查、误差判定及参数分辨率判定。
- 单次搜索的 40% 迭代预算、保持最优顶点的重初始化、按单纯形尺度估计新变化范围。
- 双线程重复搜索、单线程回退、重启时的种子偏移、两个相同最优目标值的接受条件、额外数组元素。
- 六种状态码：成功 0、初始化失败 1、中止 2、重初始化失败 3、迭代超限 4、重启超限 5。`MAX_ITERATIONS` 保留为状态 4 的兼容别名。

目标函数在定义域外应返回 NaN；这是原版优化器的算法输入，不能替换成大罚值。多线程目标函数须线程安全。原版双线程共享控制状态在 C++ 中使用原子变量和互斥锁，避免数据竞争；并行调度及不同平台数学库仍可能带来末位差异，不能保证跨语言逐位相同或任意问题全局收敛。

## 数值库的移植范围

- 不移植设置弹窗、进度显示、键盘中止或绘图界面。初值和算法配置由调用方通过数值 API 设置，
  取消拟合使用 `Minimizer::abort()`。调用方可用 `f(x)` 采样并自行绘图。
- 启用字符串公式后，`doCustomFit("y=a+b*Math.exp(c*x)", initial)` 使用 Qt `QJSEngine` 执行公式。
  支持常用数学函数，参数按原版 `a...f` 识别；公式拟合强制单线程。
  不提供 ImageJ 图像处理、窗口、插件命令或完整宏运行时。
- 报告保留原版统计项；数字文本格式使用 Qt。C++ 数组采用值语义。

数值细节也遵循上游：多项式仅解析求解偏移/斜率，其余系数参与单纯形搜索；回归 SSE 保留原版浮点误差下限；ERF 使用原版 `IJMath.erf` 近似而非系统 `erf`。`EXP_REGRESSION` / `POWER_REGRESSION` 在对数空间拟合，`RODBARD2` 交换 x/y 后拟合，因此三者的 SSE/R²/拟合优度属于变换后的空间，残差仍针对原始数据。

## 原版对照与测试

基准固定于 ImageJ 提交 [`888e13d30cf0fa7c13f601fe929fdb3f940cafea`](https://github.com/imagej/ImageJ/tree/888e13d30cf0fa7c13f601fe929fdb3f940cafea)。来源为其中的 `ij/measure/CurveFitter.java`、`ij/measure/Minimizer.java` 和 `ij/util/IJMath.java`，上游公有领域声明见 [imagej-license.txt](imagej-license.txt)。实现保留并核对了原注释中的算法理由及边界契约。

`tests/java_reference.jsonl` 包含由未修改 Java 源码生成的 86 组结果：全部 25 种模型各有默认初值、显式初值、加权三组案例，另有退化直线、负 y 指数、含原点的幂回归、NaN 清洗，以及覆盖全部六个状态的独立 Minimizer 案例。QTest 比较状态、参数、SSE、残差和统计值；六个单线程优化器案例额外精确比较迭代数与完成次数。双精度数值采用绝对/相对容差 `1e-5`，不把跨语言浮点末位视为算法偏差。

另外使用 QTest 验证并行搜索、参数分辨率、额外数组元素、自定义回调和结果报告；启用字符串公式时也验证其解析和拟合。

使用 JDK 21 可从固定提交重建并核验基准：

```sh
python3 tests/generate_reference.py --check
# 或复用固定提交的现有 ImageJ checkout：
python3 tests/generate_reference.py --source /path/to/ImageJ --check
```

脚本核对提交与三个原始文件的 SHA-256，在临时目录编译，`--check` 不写基准文件。单纯更新预期值不能代替调查算法差异。
