#pragma once

#include "curvefitterexport.h"
#include "minimizer.h"

#include <functional>
#include <QStringList>

/**
 * @brief ImageJ CurveFitter 的 Qt/C++ 移植，保留原版拟合、回归消元和统计流程。
 */
class CurveFitter
{
    Q_DECLARE_TR_FUNCTIONS(CurveFitter)

public:
    using Array = QVector<double>;
    using Parameters = Array;
    using UserFunction = std::function<double(const Array &, double)>;
    enum FitType
    {
        STRAIGHT_LINE,
        POLY2,
        POLY3,
        POLY4,
        EXPONENTIAL,
        POWER,
        LOG,
        RODBARD,
        GAMMA_VARIATE,
        LOG2,
        RODBARD2,
        EXP_WITH_OFFSET,
        GAUSSIAN,
        EXP_RECOVERY,
        INV_RODBARD,
        EXP_REGRESSION,
        POWER_REGRESSION,
        POLY5,
        POLY6,
        POLY7,
        POLY8,
        GAUSSIAN_NOOFFSET,
        EXP_RECOVERY_NOOFFSET,
        CHAPMAN,
        ERF,
        CUSTOM = 100,
    };

    enum Status
    {
        SUCCESS = Minimizer::SUCCESS,
        INITIALIZATION_FAILURE = Minimizer::INITIALIZATION_FAILURE,
        ABORTED = Minimizer::ABORTED,
        REINITIALIZATION_FAILURE = Minimizer::REINITIALIZATION_FAILURE,
        MAX_ITERATIONS_EXCEEDED = Minimizer::MAX_ITERATIONS_EXCEEDED,
        MAX_RESTARTS_EXCEEDED = Minimizer::MAX_RESTARTS_EXCEEDED,
        MAX_ITERATIONS = MAX_ITERATIONS_EXCEEDED,
    };

    CURVEFITTER_EXPORT CurveFitter(const Array &x, const Array &y);

    CURVEFITTER_EXPORT void doFit(int type);
    CURVEFITTER_EXPORT void doCustomFit(UserFunction function,
                                        int count,
                                        const QString &formula,
                                        const Array &initial = {},
                                        const Array &variations = {});
#ifdef CURVEFITTER_ENABLE_STRING_FORMULAS
    CURVEFITTER_EXPORT int doCustomFit(const QString &equation, const Array &initial = {});
#endif
    void setInitialParameters(const Array &p) { initialParams = p; }

    void setInitialParamVariations(const Array &p) { initialParamVariations = p; }

    void setWeights(const Array &w) { weights = w; }

    Minimizer &getMinimizer() { return minimizer; }

    CURVEFITTER_EXPORT void setOffsetMultiplySlopeParams(int offset, int multiply, int slope);
    CURVEFITTER_EXPORT Array getParams();
    CURVEFITTER_EXPORT Array getResiduals();
    CURVEFITTER_EXPORT double getSumResidualsSqr();
    CURVEFITTER_EXPORT double getSD();
    CURVEFITTER_EXPORT double getRSquared();
    CURVEFITTER_EXPORT double getFitGoodness();
    int getStatus() const { return linearRegressionUsed ? Minimizer::SUCCESS : minimizerStatus; }

    CURVEFITTER_EXPORT QString getStatusString() const;
    CURVEFITTER_EXPORT QString getResultString();
    void setRestarts(int n) { minimizer.setMaxRestarts(n); }

    int getRestarts() const { return minimizer.getMaxRestarts(); }

    CURVEFITTER_EXPORT void setMaxError(double error);
    int getIterations() const { return linearRegressionUsed ? 1 : minimizer.getIterations(); }

    int getMaxIterations() const { return minimizer.getMaxIterations(); }

    void setMaxIterations(int n) { minimizer.setMaxIterations(n); }

    Array getXPoints() const { return xData; }

    Array getYPoints() const { return yData; }

    int getFit() const { return fitType; }

    CURVEFITTER_EXPORT int getNumParams();
    CURVEFITTER_EXPORT QString getName();
    CURVEFITTER_EXPORT QString getFormula();
    CURVEFITTER_EXPORT QString getMacroCode();
    CURVEFITTER_EXPORT double f(double x);
    CURVEFITTER_EXPORT double f(const Array &p, double x);
    CURVEFITTER_EXPORT double userFunction(Array &p, double dummy);
    static CURVEFITTER_EXPORT double f(int type, const Array &p, double x);
    static CURVEFITTER_EXPORT int getNumParams(int type);
#ifdef CURVEFITTER_ENABLE_STRING_FORMULAS
    static CURVEFITTER_EXPORT int getNumParams(const QString &equation);
#endif
    static CURVEFITTER_EXPORT int getMax(const Array &values);
    static CURVEFITTER_EXPORT int getFitCode(const QString &name);
    static CURVEFITTER_EXPORT QStringList fitList();
    static CURVEFITTER_EXPORT QStringList formulas();
    static CURVEFITTER_EXPORT QStringList macroCodes();
    static CURVEFITTER_EXPORT QVector<int> sortedTypes();
    static CURVEFITTER_EXPORT QStringList getSortedFitList();
    static constexpr int IterFactor{500};

private:
    void minimizerParamsToFullParams(Array &p, bool regression);
    void doRegression(Array &p);
    double fullParamsToMinimizerParams(Array &p);
    void modifyInitialParamsAndVariations();
    bool makeInitialParamsAndVariations(int type);
    void getOffsetAndFactorParams();
    void calculateSumYandY2();
    bool isModifiedFitType(int type);
    bool prepareModifiedFitType(int type);
    void postProcessModifiedFitType(int type);
    static double sqr(double x) { return x * x; }

private:
    static constexpr int GAUSSIAN_INTERNAL{101}, RODBARD_INTERNAL{102};
    int fitType{-1}, numPoints{}, numParams{}, customParamCount{}, numRegressionParams{};
    int offsetParam{-1}, factorParam{-1};
    Array xData, yData, weights, xDataSave, yDataSave; // 当前拟合数据、权重及变换前的原始数据。
    double ySign{}, sumY{std::numeric_limits<double>::quiet_NaN()},
        sumY2{std::numeric_limits<double>::quiet_NaN()};
    double sumWeights{std::numeric_limits<double>::quiet_NaN()}, maxRelError{1e-10};
    qint64 time{};
    Array initialParams, initialParamVariations, minimizerInitialParams,
        minimizerInitialParamVariations, finalParams;
    QString customFormula, errorString;
    UserFunction customFunction;
    bool macroDefined{}, hasSlopeParam{}, linearRegressionUsed{}, restrictPower{};
    Minimizer minimizer;
    int minimizerStatus{Minimizer::INITIALIZATION_FAILURE};
};
