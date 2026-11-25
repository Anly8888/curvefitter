#include "curvefitter.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
////////////////////////////////////////////////////////////////////////
// class CurveFitter

namespace
{
constexpr double quietNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double inf = std::numeric_limits<double>::infinity();
int degree(int t)
{
    if (t == 0)
        return 1;
    if (t >= 1 && t <= 3)
        return t + 1;
    if (t >= 17 && t <= 20)
        return t - 12;
    return -1;
}
}
QStringList CurveFitter::fitList()
{
    return { "Straight Line", "2nd Degree Polynomial", "3rd Degree Polynomial",
        "4th Degree Polynomial", "Exponential", "Power", "Log", "Rodbard", "Gamma Variate",
        "y = a+b*ln(x-c)", "Rodbard (NIH Image)", "Exponential with Offset", "Gaussian",
        "Exponential Recovery", "Inverse Rodbard", "Exponential (linear regression)",
        "Power (linear regression)", "5th Degree Polynomial", "6th Degree Polynomial",
        "7th Degree Polynomial", "8th Degree Polynomial", "Gaussian (no offset)",
        "Exponential Recovery (no offset)", "Chapman-Richards", "Error Function" };
}
QStringList CurveFitter::formulas()
{
    return { "y = a+bx", "y = a+bx+cx^2", "y = a+bx+cx^2+dx^3", "y = a+bx+cx^2+dx^3+ex^4",
        "y = a*exp(bx)", "y = a*x^b", "y = a*ln(bx)", "y = d+(a-d)/(1+(x/c)^b)",
        "y = b*(x-a)^c*exp(-(x-a)/d)", "y = a+b*ln(x-c)",
        "x = d+(a-d)/(1+(y/c)^b) [y = c*((x-a)/(d-x))^(1/b)]", "y = a*exp(-bx) + c",
        "y = a + (b-a)*exp(-(x-c)*(x-c)/(2*d*d))", "y = a*(1-exp(-b*x)) + c",
        "y = c*((x-a)/(d-x))^(1/b)", "y = a*exp(bx)", "y = a*x^b", "y = a+bx+cx^2+dx^3+ex^4+fx^5",
        "y = a+bx+cx^2+dx^3+ex^4+fx^5+gx^6", "y = a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7",
        "y = a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7+ix^8", "y = a*exp(-(x-b)*(x-b)/(2*c*c))",
        "y = a*(1-exp(-b*x))", "y = a*(1-exp(-b*x))^c", "y = a+b*erf((x-c)/d)" };
}
QVector<int> CurveFitter::sortedTypes()
{
    return { 0, 1, 2, 3, 17, 18, 19, 20, 5, 16, 4, 15, 11, 13, 22, 6, 9, 12, 21, 24, 7, 10, 14, 8,
        23 };
}
QStringList CurveFitter::getSortedFitList()
{
    QStringList r;
    for (int t : sortedTypes())
        r << fitList()[t];
    return r;
}
int CurveFitter::getFitCode(const QString& s)
{
    return fitList().indexOf(s);
}
int CurveFitter::getNumParams(int t)
{
    if (degree(t) >= 0)
        return degree(t) + 1;
    switch (t)
    {
    case EXPONENTIAL:
    case POWER:
    case LOG:
    case EXP_REGRESSION:
    case POWER_REGRESSION:
    case EXP_RECOVERY_NOOFFSET:
        return 2;
    case LOG2:
    case EXP_WITH_OFFSET:
    case EXP_RECOVERY:
    case GAUSSIAN_NOOFFSET:
    case CHAPMAN:
        return 3;
    case RODBARD:
    case RODBARD2:
    case INV_RODBARD:
    case GAMMA_VARIATE:
    case GAUSSIAN:
    case ERF:
        return 4;
    default:
        throw std::invalid_argument("Invalid fit type");
    }
}
int CurveFitter::getNumParams() const
{
    return type_ == CUSTOM ? customCount_ : getNumParams(type_);
}
QString CurveFitter::getName() const
{
    return type_ == CUSTOM ? "Custom" : type_ >= 0 ? fitList().value(type_) : QString();
}
QString CurveFitter::getFormula() const
{
    return type_ == CUSTOM ? customFormula_ : formulas().value(type_);
}
double CurveFitter::f(const Parameters& p, double x) const
{
    return type_ == CUSTOM ? custom_(p, x) : f(type_, p, x);
}
double CurveFitter::f(int t, const Parameters& p, double x)
{
    if (p.size() < getNumParams(t))
        throw std::invalid_argument("Too few parameters");
    if (degree(t) >= 0)
    {
        double y = p[degree(t)];
        for (int i = degree(t) - 1; i >= 0; --i)
            y = y * x + p[i];
        return y;
    }
    switch (t)
    {
    case EXPONENTIAL:
    case EXP_REGRESSION:
        return p[0] * std::exp(p[1] * x);
    case POWER:
    case POWER_REGRESSION:
        return p[0] * std::pow(x, p[1]);
    case LOG:
        return x == 0 ? -1000 * p[0] : p[0] * std::log(p[1] * x);
    case EXP_WITH_OFFSET:
        return p[0] * std::exp(-p[1] * x) + p[2];
    case EXP_RECOVERY:
        return p[0] * (-std::expm1(-p[1] * x)) + p[2];
    case EXP_RECOVERY_NOOFFSET:
        return p[0] * (-std::expm1(-p[1] * x));
    case CHAPMAN:
        return p[0] * std::pow(-std::expm1(-p[1] * x), p[2]);
    case GAUSSIAN:
        return p[0] + (p[1] - p[0]) * std::exp(-std::pow((x - p[2]) / p[3], 2) / 2);
    case GAUSSIAN_NOOFFSET:
        return p[0] * std::exp(-std::pow((x - p[1]) / p[2], 2) / 2);
    case RODBARD:
        return p[3] + (p[0] - p[3]) / (1 + std::pow(x / p[2], p[1]));
    case GAMMA_VARIATE:
        if (p[0] >= x)
            return 0;
        if (p[1] <= 0 || p[2] <= 0 || p[3] <= 0)
            return quietNaN;
        return p[1] * std::pow(x - p[0], p[2]) * std::exp(-(x - p[0]) / p[3]);
    case LOG2:
        return x <= p[2] ? quietNaN : p[0] + p[1] * std::log(x - p[2]);
    case INV_RODBARD:
    case RODBARD2:
        if (p[3] - x < 2 * std::numeric_limits<double>::denorm_min() || x < p[0])
            return t == RODBARD2 ? 0 : quietNaN;
        return p[2] * std::pow((x - p[0]) / (p[3] - x), 1 / p[1]);
    case ERF:
        return p[0] + p[1] * std::erf((x - p[2]) / p[3]);
    default:
        return quietNaN;
    }
}
CurveFitter::CurveFitter(const QVector<double>& x, const QVector<double>& y)
    : originalCount_(x.size())
{
    if (x.size() != y.size())
        throw std::invalid_argument("x/y size mismatch");
    for (int i = 0; i < x.size(); ++i)
        if (!std::isnan(x[i] + y[i]))
        {
            x_ << x[i];
            y_ << y[i];
            kept_ << i;
        }
    weights_.fill(1, x_.size());
}
/**
 * @brief 设置原始输入各点的权重；零权重点不参与拟合。
 * @param w 非负有限权重，长度须等于过滤 NaN 前的输入长度。
 */
void CurveFitter::setWeights(const QVector<double>& w)
{
    if (w.size() != originalCount_)
        throw std::invalid_argument("Weights size mismatch");
    for (double v : w)
        if (!std::isfinite(v) || v < 0)
            throw std::invalid_argument("Weights must be finite and nonnegative");
    for (int i = 0; i < kept_.size(); ++i)
        weights_[i] = w[kept_[i]];
}
void CurveFitter::setMaxIterations(int n)
{
    if (n < 1)
        throw std::invalid_argument("Iterations must be positive");
    maxIterations_ = n;
}
void CurveFitter::setRestarts(int n)
{
    if (n < 0)
        throw std::invalid_argument("Restarts must be nonnegative");
    restarts_ = n;
}
void CurveFitter::setMaxError(double e)
{
    if (!std::isfinite(e) || e <= 0)
        throw std::invalid_argument("Tolerance must be positive");
    tolerance_ = e;
}
/**
 * @brief 使用内置模型拟合；通过状态接口区分成功、无效数据及迭代上限。
 * @param t ImageJ 模型编号。
 */
void CurveFitter::doFit(int t)
{
    getNumParams(t);
    type_ = t;
    fit();
}
/**
 * @brief 用 C++ 回调替代 ImageJ UserFunction；初值必须覆盖所有参数。
 * @param fn 模型函数，异常直接传给调用者。
 * @param n 参数个数。
 * @param formula 用于结果报告的公式文本。
 * @param p 初始参数。
 * @param v 初始搜索步长，省略时按参数尺度估计。
 */
void CurveFitter::doCustomFit(
    UserFunction fn, int n, const QString& formula, const Parameters& p, const Parameters& v)
{
    if (!fn || n < 1 || p.size() != n)
        throw std::invalid_argument("Invalid custom function");
    custom_ = std::move(fn);
    customCount_ = n;
    customFormula_ = formula;
    initial_ = p;
    variations_ = v;
    type_ = CUSTOM;
    fit();
}

namespace
{
// Householder QR avoids squaring the condition number as normal equations do.
bool leastSquares(QVector<QVector<double>> a, QVector<double> b, QVector<double>& solution)
{
    const int rows = a.size(), cols = a[0].size();
    if (rows < cols)
        return false;
    for (int k = 0; k < cols; ++k)
    {
        double norm = 0;
        for (int i = k; i < rows; ++i)
            norm = std::hypot(norm, a[i][k]);
        if (norm == 0 || !std::isfinite(norm))
            return false;
        double alpha = -std::copysign(norm, a[k][k]);
        QVector<double> v(rows - k);
        for (int i = k; i < rows; ++i)
            v[i - k] = a[i][k];
        v[0] -= alpha;
        double vv = 0;
        for (double z : v)
            vv += z * z;
        for (int j = k; j < cols; ++j)
        {
            double dot = 0;
            for (int i = k; i < rows; ++i)
                dot += v[i - k] * a[i][j];
            for (int i = k; i < rows; ++i)
                a[i][j] -= 2 * v[i - k] * dot / vv;
        }
        double dot = 0;
        for (int i = k; i < rows; ++i)
            dot += v[i - k] * b[i];
        for (int i = k; i < rows; ++i)
            b[i] -= 2 * v[i - k] * dot / vv;
    }
    solution.fill(0, cols);
    for (int i = cols - 1; i >= 0; --i)
    {
        double z = b[i];
        for (int j = i + 1; j < cols; ++j)
            z -= a[i][j] * solution[j];
        solution[i] = z / a[i][i];
        if (!std::isfinite(solution[i]))
            return false;
    }
    return true;
}
}

void CurveFitter::fit()
{
    status_ = INITIALIZATION_FAILURE;
    message_ = "Invalid data or initial parameters";
    iterations_ = 0;
    const int n = getNumParams();
    params_.fill(quietNaN, n + 1);
    sse_ = quietNaN;
    variance_ = 0;
    objectiveCount_ = 0;
    QVector<double> x, y, w;
    bool restrictPower = false;
    double sign = 0;
    for (int i = 0; i < x_.size(); ++i)
    {
        if (!std::isfinite(x_[i]) || !std::isfinite(y_[i]))
            return;
        if (weights_[i] == 0)
            continue;
        double xx = x_[i], yy = y_[i];
        if (type_ == POWER_REGRESSION || type_ == EXP_REGRESSION)
        {
            if (type_ == POWER_REGRESSION)
            {
                if (xx == 0 && yy == 0)
                {
                    restrictPower = true;
                    continue;
                }
                if (xx <= 0)
                    return;
                xx = std::log(xx);
            }
            if (sign == 0)
                sign = std::copysign(1., yy);
            if (yy * sign <= 0)
                return;
            yy = std::log(yy * sign);
        }
        else if (type_ == RODBARD2)
        {
            std::swap(xx, yy);
        }
        x << xx;
        y << yy;
        w << weights_[i];
    }
    if (x.size() < n)
        return;
    double sw = 0, mean = 0;
    for (int i = 0; i < y.size(); ++i)
    {
        sw += w[i];
        mean += w[i] * y[i];
    }
    mean /= sw;
    for (int i = 0; i < y.size(); ++i)
        variance_ += w[i] * (y[i] - mean) * (y[i] - mean);
    objectiveCount_ = y.size();
    const int model = type_ == RODBARD2 ? RODBARD : type_;
    const bool logRegression = type_ == POWER_REGRESSION || type_ == EXP_REGRESSION;
    const int deg = logRegression ? 1 : degree(model);
    if (deg >= 0)
    {
        // Center and scale before QR; transform coefficients back to the public x basis.
        double center
            = (*std::min_element(x.begin(), x.end()) + *std::max_element(x.begin(), x.end())) / 2;
        double scale = 0;
        for (double v : x)
            scale = std::max(scale, std::abs(v - center));
        if (scale == 0)
            return;
        QVector<QVector<double>> a;
        QVector<double> b, c;
        for (int i = 0; i < x.size(); ++i)
        {
            QVector<double> row;
            double v = std::sqrt(w[i]);
            for (int j = 0; j <= deg; ++j)
            {
                row << v;
                v *= (x[i] - center) / scale;
            }
            a << row;
            b << std::sqrt(w[i]) * y[i];
        }
        if (!leastSquares(a, b, c))
            return;
        params_.fill(0, n + 1);
        QVector<double> power { 1 };
        for (int j = 0; j <= deg; ++j)
        {
            for (int k = 0; k < power.size(); ++k)
                params_[k] += c[j] * power[k];
            QVector<double> next(power.size() + 1, 0);
            for (int k = 0; k < power.size(); ++k)
            {
                next[k] -= power[k] * center / scale;
                next[k + 1] += power[k] / scale;
            }
            power = next;
        }
        if (restrictPower && params_[1] <= 0)
        {
            message_ = "Power must be positive when (0,0) is present";
            params_.fill(quietNaN, n + 1);
            return;
        }
        sse_ = 0;
        for (int i = 0; i < x.size(); ++i)
        {
            double pred = 0;
            for (int j = deg; j >= 0; --j)
                pred = pred * x[i] + params_[j];
            sse_ += w[i] * std::pow(pred - y[i], 2);
        }
        if (logRegression)
            params_[0] = sign * std::exp(params_[0]);
        params_[n] = sse_;
        status_ = SUCCESS;
        message_ = "Success";
        return;
    }
    double xmin = *std::min_element(x.begin(), x.end()),
           xmax = *std::max_element(x.begin(), x.end());
    double ymin = *std::min_element(y.begin(), y.end()),
           ymax = *std::max_element(y.begin(), y.end());
    double range = xmax - xmin;
    if (range <= 0)
        return;
    if ((model == POWER || model == CHAPMAN) && xmin < 0)
        return;
    if (model == LOG && xmin <= 0)
        return;
    if (model == RODBARD && xmin < 0 && xmax > 0)
        return;
    int offset = -1, factor = -1;
    switch (model)
    {
    case EXPONENTIAL:
    case POWER:
    case LOG:
    case CHAPMAN:
    case EXP_RECOVERY_NOOFFSET:
    case GAUSSIAN_NOOFFSET:
        factor = 0;
        break;
    case EXP_WITH_OFFSET:
    case EXP_RECOVERY:
        offset = 2;
        factor = 0;
        break;
    case LOG2:
    case GAUSSIAN:
    case ERF:
        offset = 0;
        factor = 1;
        break;
    case RODBARD:
        offset = 3;
        factor = 0;
        break;
    case INV_RODBARD:
        factor = 2;
        break;
    case GAMMA_VARIATE:
        factor = 1;
        break;
    }
    Parameters p(n, 1);
    const double peak = x[std::distance(y.begin(), std::max_element(y.begin(), y.end()))];
    switch (model)
    {
    case EXPONENTIAL:
        p = { mean, (y.last() > y.first() ? 1. : -1.) / range };
        break;
    case EXP_WITH_OFFSET:
    case EXP_RECOVERY:
        p = { ymax - ymin, 1 / range, ymin };
        break;
    case EXP_RECOVERY_NOOFFSET:
        p = { ymax, 1 / range };
        break;
    case POWER:
        p = { 1, 1 };
        break;
    case LOG:
        p = { mean, std::exp(1.) / xmax };
        break;
    case LOG2:
        p = { mean, 1, xmin - .1 * range };
        break;
    case RODBARD:
        p = { y.first(), 1, xmin < 0 ? xmin : xmax, y.last() };
        break;
    case INV_RODBARD:
        p = { xmin - .1 * range, 1, ymax, xmax + range };
        break;
    case GAMMA_VARIATE:
        p = { xmin, 1, std::sqrt(std::max(peak - xmin, .1 * range)),
            std::sqrt(std::max(peak - xmin, .1 * range)) };
        break;
    case GAUSSIAN:
        p = { ymin, ymax, peak,
            std::max(.01 * range, .39894 * range * (mean - ymin) / (ymax - ymin + 1e-100)) };
        break;
    case GAUSSIAN_NOOFFSET:
        p = { ymax, peak, .2 * range };
        break;
    case CHAPMAN:
        p = { ymax, 2 / range, 1.5 };
        break;
    case ERF:
        p = { mean, (ymax - ymin) / 2, (xmin + xmax) / 2, .1 * range };
        break;
    }
    if (!initial_.isEmpty())
    {
        if (initial_.size() != n)
            throw std::invalid_argument("Initial parameter count mismatch");
        p = initial_;
    }
    if (!variations_.isEmpty() && variations_.size() != n)
        throw std::invalid_argument("Variation count mismatch");
    for (double v : p)
        if (!std::isfinite(v))
            return;
    for (double v : variations_)
        if (!std::isfinite(v) || v <= 0)
            throw std::invalid_argument("Variations must be positive");
    // Internally separate the additive baseline and multiplicative amplitude.
    if (model == GAUSSIAN)
        p[1] -= p[0];
    if (model == RODBARD)
        p[0] -= p[3];
    QVector<int> free;
    for (int j = 0; j < n; ++j)
        if (j != offset && j != factor)
            free << j;
    auto value = [&](Parameters q, double xx)
    {
        if (model == CUSTOM)
            return custom_(q, xx);
        if (model == GAUSSIAN)
            q[1] += q[0];
        if (model == RODBARD)
            q[0] += q[3];
        return CurveFitter::f(model, q, xx);
    };
    auto objective = [&](Parameters& q)
    {
        if (factor >= 0)
        {
            Parameters basis = q;
            basis[factor] = 1;
            if (offset >= 0)
                basis[offset] = 0;
            double sg = 0, sgg = 0, sgy = 0, sy = 0;
            for (int i = 0; i < x.size(); ++i)
            {
                double z = value(basis, x[i]);
                if (!std::isfinite(z))
                    return inf;
                sg += w[i] * z;
                sgg += w[i] * z * z;
                sgy += w[i] * z * y[i];
                sy += w[i] * y[i];
            }
            double denom = offset >= 0 ? sgg - sg * sg / sw : sgg;
            if (denom <= 0)
                return inf;
            q[factor] = (offset >= 0 ? sgy - sg * sy / sw : sgy) / denom;
            if (offset >= 0)
                q[offset] = (sy - q[factor] * sg) / sw;
        }
        double cost = 0;
        for (int i = 0; i < x.size(); ++i)
        {
            double z = value(q, x[i]);
            if (!std::isfinite(z))
                return inf;
            cost += w[i] * (z - y[i]) * (z - y[i]);
        }
        return std::isfinite(cost) ? cost : inf;
    };
    QVector<double> steps;
    for (int j : free)
        steps << (!variations_.isEmpty() ? variations_[j]
                                         : std::max(std::abs(p[j]) * .1, range * .01));
    Parameters best = p;
    double bestCost = objective(best);
    bool converged = false;
    const int dim = free.size();
    for (int restart = 0; restart <= restarts_ && iterations_ < maxIterations_; ++restart)
    {
        QVector<Parameters> vertices(dim + 1, best);
        QVector<double> costs(dim + 1);
        for (int k = 1; k <= dim; ++k)
            vertices[k][free[k - 1]] += steps[k - 1] * std::pow(.5, restart);
        for (int k = 0; k <= dim; ++k)
            costs[k] = objective(vertices[k]);
        while (iterations_ < maxIterations_)
        {
            ++iterations_;
            QVector<int> order;
            for (int k = 0; k <= dim; ++k)
                order << k;
            std::sort(
                order.begin(), order.end(), [&](int a, int b) { return costs[a] < costs[b]; });
            int lo = order[0], hi = order[dim], second = order[dim - 1];
            if (costs[lo] < bestCost)
            {
                bestCost = costs[lo];
                best = vertices[lo];
            }
            double spread = 0;
            for (int j : free)
                spread = std::max(spread,
                    std::abs(vertices[hi][j] - vertices[lo][j]) / (1 + std::abs(vertices[lo][j])));
            if (std::isfinite(costs[lo])
                && costs[hi] - costs[lo] <= tolerance_ * std::max(1., costs[lo])
                && spread < std::sqrt(tolerance_))
            {
                converged = true;
                break;
            }
            Parameters center = vertices[lo];
            for (int j : free)
            {
                center[j] = 0;
                for (int k = 0; k <= dim; ++k)
                    if (k != hi)
                        center[j] += vertices[k][j] / dim;
            }
            auto trial = [&](double factor)
            {
                Parameters q = center;
                for (int j : free)
                    q[j] += factor * (center[j] - vertices[hi][j]);
                return q;
            };
            Parameters r = trial(1);
            double cr = objective(r);
            if (cr < costs[lo])
            {
                Parameters e = trial(2);
                double ce = objective(e);
                vertices[hi] = ce < cr ? e : r;
                costs[hi] = std::min(ce, cr);
            }
            else if (cr < costs[second])
            {
                vertices[hi] = r;
                costs[hi] = cr;
            }
            else
            {
                Parameters c = trial(cr < costs[hi] ? .5 : -.5);
                double cc = objective(c);
                if (cc < std::min(cr, costs[hi]))
                {
                    vertices[hi] = c;
                    costs[hi] = cc;
                }
                else
                    for (int k = 0; k <= dim; ++k)
                        if (k != lo)
                        {
                            for (int j : free)
                                vertices[k][j] = (vertices[k][j] + vertices[lo][j]) / 2;
                            costs[k] = objective(vertices[k]);
                        }
            }
        }

        // 达到迭代上限时也保留最后一次试探找到的更优解。
        for (int k = 0; k <= dim; ++k)
        {
            if (costs[k] < bestCost)
            {
                bestCost = costs[k];
                best = vertices[k];
            }
        }
    }

    for (double v : best)
        if (!std::isfinite(v))
            return;
    if (!std::isfinite(bestCost))
        return;
    p = best;
    if (model == GAUSSIAN)
    {
        p[1] += p[0];
        p[3] = std::abs(p[3]);
    }
    if (model == RODBARD)
        p[0] += p[3];
    if (model == GAUSSIAN_NOOFFSET)
        p[2] = std::abs(p[2]);
    if (model == ERF && p[3] < 0)
    {
        p[1] = -p[1];
        p[3] = -p[3];
    }
    sse_ = bestCost;
    params_ = p;
    params_ << sse_;
    status_ = converged ? SUCCESS : MAX_ITERATIONS;
    message_ = converged ? "Success" : "Maximum iterations reached";
}

QVector<double> CurveFitter::getResiduals() const
{
    QVector<double> r;
    for (int i = 0; i < x_.size(); ++i)
        r << y_[i] - f(params_, x_[i]);
    return r;
}
double CurveFitter::getRSquared() const
{
    return std::isnan(sse_) ? quietNaN : variance_ > 0 ? 1 - sse_ / variance_ : 0;
}
double CurveFitter::getFitGoodness() const
{
    int df = objectiveCount_ - getNumParams();
    return std::isnan(sse_)       ? quietNaN
        : variance_ > 0 && df > 0 ? 1 - sse_ / variance_ * objectiveCount_ / df
                                  : 0;
}
double CurveFitter::getSD() const
{
    auto r = getResiduals();
    if (r.size() < 2)
        return quietNaN;
    double mean = std::accumulate(r.begin(), r.end(), 0.) / r.size(), v = 0;
    for (double z : r)
        v += (z - mean) * (z - mean);
    return std::sqrt(v / (r.size() - 1));
}
QString CurveFitter::getResultString() const
{
    QString r = getName() + "\n" + getFormula() + "\n" + message_;
    for (int i = 0; i < getNumParams(); ++i)
        r += QString("\n%1 = %2").arg(QChar('a' + i)).arg(params_.value(i, quietNaN), 0, 'g', 12);
    return r + QString("\nSSE = %1\nR² = %2").arg(sse_, 0, 'g', 12).arg(getRSquared(), 0, 'g', 12);
}
