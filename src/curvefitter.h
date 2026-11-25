#pragma once
#include <QStringList>
#include <QVector>
#include <functional>
#include <limits>

/**
 * @brief ImageJ 拟合模型的 Qt Core 数值库，支持加权回归和自定义函数。
 */
class CurveFitter
{
public:
    /**
     * @brief 保留 ImageJ 内置模型编号及参数顺序。
     */
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

    /**
     * @brief 本库的拟合状态，数值不对应 Java Minimizer 状态。
     */
    enum Status
    {
        SUCCESS,
        INITIALIZATION_FAILURE,
        MAX_ITERATIONS,
    };

    using Parameters = QVector<double>;
    using UserFunction = std::function<double(const Parameters&, double)>;

    CurveFitter(const QVector<double>& x, const QVector<double>& y);

    void doFit(int type);
    void doCustomFit(UserFunction function, int count, const QString& formula,
        const Parameters& initial, const Parameters& variations = { });

    void setInitialParameters(const Parameters& p) { initial_ = p; }

    void setInitialParamVariations(const Parameters& p) { variations_ = p; }

    // Weights refer to the original input, including any NaN pairs removed by the constructor.
    void setWeights(const QVector<double>& weights);
    void setMaxIterations(int n);
    void setMaxError(double tolerance);
    void setRestarts(int n);

    Parameters getParams() const { return params_; } // parameters followed by objective SSE

    QVector<double> getResiduals() const;

    QVector<double> getXPoints() const { return x_; }

    QVector<double> getYPoints() const { return y_; }

    double getSumResidualsSqr() const { return sse_; }

    double getRSquared() const;
    double getFitGoodness() const;
    double getSD() const;

    int getStatus() const { return status_; }

    QString getStatusString() const { return message_; }

    int getIterations() const { return iterations_; }

    int getNumParams() const;
    QString getName() const;
    QString getFormula() const;
    QString getResultString() const;
    double f(const Parameters& p, double x) const;

    static double f(int type, const Parameters& p, double x);
    static int getNumParams(int type);
    static int getFitCode(const QString& name);
    static QStringList fitList();
    static QStringList formulas();
    static QVector<int> sortedTypes();
    static QStringList getSortedFitList();

private:
    void fit();

private:
    QVector<double> x_, y_, weights_;
    QVector<int> kept_;
    int originalCount_ { }, type_ { -1 }, customCount_ { }, status_ { INITIALIZATION_FAILURE };
    int iterations_ { }, maxIterations_ { 20000 }, restarts_ { 2 }, objectiveCount_ { };
    double tolerance_ { 1e-10 }, sse_ { std::numeric_limits<double>::quiet_NaN() }, variance_ { };
    Parameters initial_, variations_, params_;
    UserFunction custom_;
    QString customFormula_, message_ { "No fit performed" };
};
