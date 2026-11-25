#pragma once

#include <atomic>
#include <cstdint>
#include <functional>
#include <mutex>
#include <QCoreApplication>
#include <QElapsedTimer>
#include <QThread>
#include <QVector>

/**
 * @brief ImageJ Minimizer 的随机重启 Nelder–Mead 实现，保留上游状态码和搜索流程。
 *
 * 保留 Lagarias 等人的外收缩扩展、最佳顶点中心检查，以及退化单纯形重初始化。
 * 重复搜索在两个最优目标值满足误差限制后结束，但不能保证所得最小值是全局最小值。
 * 目标函数的有效定义域应为凸集；定义域外返回 NaN 而非罚值有助于可靠初始化。
 */
class Minimizer
{
    Q_DECLARE_TR_FUNCTIONS(Minimizer)

public:
    using Array = QVector<double>;
    using Matrix = QVector<Array>;
    using UserFunction = std::function<double(Array &, double)>;
    using StatusCallback = std::function<void(const QString &)>;
    using EscapeCallback = std::function<bool()>;

    enum Status
    {
        SUCCESS,                  // 成功收敛。
        INITIALIZATION_FAILURE,   // 无法找到有效初始单纯形。
        ABORTED,                  // 调用方或 ESC 请求中止。
        REINITIALIZATION_FAILURE, // 无法重新构造有效单纯形，结果可能不准确。
        MAX_ITERATIONS_EXCEEDED,  // 达到搜索迭代预算。
        MAX_RESTARTS_EXCEEDED,    // 重启预算用完，未找到两个一致的最优值。
    };

    Minimizer() = default;

    void setFunction(UserFunction function, int count);
    int minimize(const Array &initialParams = {}, const Array &variations = {});
    int minimizeOnce(const Array &initialParams = {}, const Array &variations = {});
    Array getParams();
    double getFunctionValue();
    int getIterations() const { return totalNumIter; }

    void setMaxIterations(int n) { maxIter = n; }

    int getMaxIterations() const { return maxIter; }

    void setMaxRestarts(int n) { maxRestarts = n; }

    int getMaxRestarts() const { return maxRestarts; }

    int getCompletedMinimizations() const { return numCompletedMinimizations; }

    void setRandomSeed(int seed) { randomSeed = seed; }

    void setMaxError(double relative) { maxRelError = relative; }

    void setMaxError(double relative, double absolute);
    void setParamResolutions(const Array &resolutions) { paramResolutions = resolutions; }

    void setMaximumThreads(int n) { useSingleThread = n <= 1; }

    void abort() { status = ABORTED; }

    void setStatusAndEsc(const QString &prefix, bool escape);
    void setStatusCallbacks(StatusCallback display, EscapeCallback escape = {});
    void setExtraArrayElements(int count) { numExtraArrayElements = count; }

    static QStringList statusStrings();

private:
    /**
     * @brief 与 java.util.Random 相同的 48 位线性同余发生器。
     */
    class JavaRandom
    {
    public:
        explicit JavaRandom(int seed);
        float nextFloat();
        double nextDouble();

    private:
        std::uint32_t next(int bits);
        std::uint64_t seed{};
    };

    void minimizeOnce(const Array &initialParams, const Array &variations, int seed);
    void reportProgress();
    static void copyVertex(const Array &source, Array &target) { target = source; }

    int minimizeSimplex(Matrix &simp);
    void getVertexAndEvaluate(Array &center, Array &worstVertex, double howFar, Array &newVertex);
    void getCenter(Matrix &simp, int excludeVertex, Array &center);
    void shrinkSimplexAndEvaluate(Matrix &simp, int best);
    bool belowErrorLimit(double highest, double lowest, double sensitivity);
    Matrix makeSimplex(const Array &initialParams,
                       const Array &initialParamVariations,
                       JavaRandom &random);
    bool belowResolutionLimit(Matrix &simp, int best);
    bool belowResolutionLimit(Array &vertex1, Array &vertex2);
    void findValidInitalParams(Array &params,
                               const Array &initialParamVariations,
                               JavaRandom &random);
    bool reInitializeSimplex(Matrix &simp,
                             int bestVertexNumber,
                             const Array &paramVariations,
                             JavaRandom &random);
    bool initializeSimplex(Matrix &simp, const Array &paramVariations, JavaRandom &random);
    Array makeNewParamVariations(Matrix &simp,
                                 int bestVertexNumber,
                                 const Array &initialParams,
                                 const Array &initialParamVariations);
    void evaluate(Array &vertex);
    double value(Array &vertex);
    void order(Matrix &simp, QVector<int> &worstNextBestArray);

private:
    static constexpr double C_REFLECTION{1.0}, C_CONTRACTION{0.5}, C_EXPANSION{2.0}, C_SHRINK{0.5};
    static constexpr int ITER_FACTOR{750}, WORST{}, NEXT_WORST{1}, BEST{2};
    int numParams{}, numVertices{}, numExtraArrayElements{};
    UserFunction userFunction; // 允许通过额外数组元素传回辅助结果的目标函数。
    double maxRelError{1e-10}, maxAbsError{1e-100};
    Array paramResolutions;    // 各参数可分辨的最小变化，空数组表示不检查。
    int maxIter{}, maxRestarts{2}, randomSeed{};
    std::atomic<int> totalNumIter{}, numCompletedMinimizations{}, status{SUCCESS};
    std::atomic<bool> wasInitialized{};
    bool useSingleThread{QThread::idealThreadCount() <= 1};
    Array result;         // 最优参数、目标值及用户保留的额外元素。
    Matrix resultsVector; // 多次搜索的结果，主线程在 join 后筛选。
    bool collectResults{};
    QString ijStatusString;
    bool checkEscape{};
    int nextIterationForStatus{10};
    QElapsedTimer timer;
    StatusCallback statusCallback;
    EscapeCallback escapeCallback;
    std::mutex resultsMutex, progressMutex; // 分别保护结果提交和宿主状态回调。
};
