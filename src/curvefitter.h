#pragma once

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

    CurveFitter(const Array &x, const Array &y);

    void doFit(int type);
    void doCustomFit(UserFunction function,
                     int count,
                     const QString &formula,
                     const Array &initial = {},
                     const Array &variations = {});
#ifdef CURVEFITTER_ENABLE_STRING_FORMULAS
    int doCustomFit(const QString &equation, const Array &initial = {});
#endif
    void setInitialParameters(const Array &p) { initialParams = p; }

    void setInitialParamVariations(const Array &p) { initialParamVariations = p; }

    void setWeights(const Array &w) { weights = w; }

    Minimizer &getMinimizer() { return minimizer; }

    void setOffsetMultiplySlopeParams(int offset, int multiply, int slope);
    Array getParams();
    Array getResiduals();
    double getSumResidualsSqr();
    double getSD();
    double getRSquared();
    double getFitGoodness();
    int getStatus() const { return linearRegressionUsed ? Minimizer::SUCCESS : minimizerStatus; }

    QString getStatusString() const;
    QString getResultString();
    void setRestarts(int n) { minimizer.setMaxRestarts(n); }

    int getRestarts() const { return minimizer.getMaxRestarts(); }

    void setMaxError(double error);
    int getIterations() const { return linearRegressionUsed ? 1 : minimizer.getIterations(); }

    int getMaxIterations() const { return minimizer.getMaxIterations(); }

    void setMaxIterations(int n) { minimizer.setMaxIterations(n); }

    Array getXPoints() const { return xData; }

    Array getYPoints() const { return yData; }

    int getFit() const { return fitType; }

    int getNumParams();
    QString getName();
    QString getFormula();
    QString getMacroCode();
    double f(double x);
    double f(const Array &p, double x);
    double userFunction(Array &p, double dummy);
    static double f(int type, const Array &p, double x);
    static int getNumParams(int type);
#ifdef CURVEFITTER_ENABLE_STRING_FORMULAS
    static int getNumParams(const QString &equation);
#endif
    static int getMax(const Array &values);
    static int getFitCode(const QString &name);
    static QStringList fitList();
    static QStringList formulas();
    static QStringList macroCodes();
    static QVector<int> sortedTypes();
    static QStringList getSortedFitList();
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
