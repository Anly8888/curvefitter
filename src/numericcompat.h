#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <QVector>

namespace ImageJNumeric {

inline constexpr double qNaN = std::numeric_limits<double>::quiet_NaN();

template<class A, class B>
auto jmin(A a, B b)
{
    using T = std::common_type_t<A, B>;
    if (std::isnan(static_cast<double>(a)) || std::isnan(static_cast<double>(b)))
        return static_cast<T>(qNaN);
    return std::min<T>(a, b);
}

template<class A, class B>
auto jmax(A a, B b)
{
    using T = std::common_type_t<A, B>;
    if (std::isnan(static_cast<double>(a)) || std::isnan(static_cast<double>(b)))
        return static_cast<T>(qNaN);
    return std::max<T>(a, b);
}

inline double signum(double x)
{
    return std::isnan(x) || x == 0 ? x : x > 0 ? 1 : -1;
}

inline void fill(QVector<double> &a, double x)
{
    std::fill(a.begin(), a.end(), x);
}

inline void arraycopy(
    const QVector<double> &source, int start, QVector<double> &target, int dest, int count)
{
    // Java arraycopy permits overlapping ranges and copies values before overwriting.
    const auto copy = source.mid(start, count);
    for (int i = 0; i < count; ++i)
        target[dest + i] = copy[i];
}

/**
 * @brief 与 ImageJ IJMath.erf 相同的近似式，避免替换为系统 erf 改变搜索路径。
 * @param x 自变量。
 * @return 最大绝对和相对误差小于 1e-12 的误差函数近似。
 */
inline double imagejErf(double x)
{
    constexpr double A = .05857895078654250866288, B = -.00626245895772819579,
                     C = -.00299946450696036814, D = .289389696496082416, E = .0539962589851632982,
                     F = .00508516909930653109, G = .000215969713046142876,
                     H = -.000225663858340491571, I = -3.06833213472529049e-7;
    const double x2 = x * x, twoOverSqrtPi = 2 / std::sqrt(3.14159265358979323846);
    if (x2 < 1e-8)
        return twoOverSqrtPi * x * (1 + x2 * (-1. / 3. + x2 * (1. / 10.)));
    const double v = x
                     * ((twoOverSqrtPi - A)
                        + A * (1 + x2 * x2 * (B + x2 * (C + x2 * (H + x2 * I))))
                              / (1 + x2 * (D + x2 * (E + x2 * (F + x2 * G)))));
    const double erf = x2 > 36 ? 1 : std::sqrt(1 - std::exp(-v * v));
    return x > 0 ? erf : -erf;
}

} // namespace ImageJNumeric
