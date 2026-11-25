#include "curvefitter.h"
#include "numericcompat.h"

#include <memory>
#include <QJSEngine>
#include <QJSValue>
#include <QRegularExpression>

////////////////////////////////////////////////////////////////////////
// class CurveFitter

/**
 * @brief 按 ImageJ 的 a...f 约定统计自定义公式参数。
 * @param equation 必须包含 x 与 y 的公式文本。
 * @return 出现的参数个数，缺少 x 或 y 时返回 0。
 */
int CurveFitter::getNumParams(const QString &equation)
{
    QString code = equation;
    code.remove(QRegularExpression("/\\*[\\s\\S]*?\\*/|//[^\\n]*|\"(?:\\\\.|[^\"\\\\])*\""));
    const auto contains = [&](const QString &word) {
        return QRegularExpression("\\b" + word + "\\b").match(code).hasMatch();
    };
    if (!contains("x") || !contains("y"))
        return 0;
    int count{};
    for (char parameter = 'a'; parameter <= 'f'; ++parameter)
        if (contains(QString(QChar(parameter))))
            ++count;

    return count;
}

/**
 * @brief 使用 Qt JavaScript 引擎执行数值公式，替代 ImageJ 宏解释器的宿主接口。
 * @param equation y=... 形式的数值公式；支持 Math 函数和常用裸函数名。
 * @param initial 初始参数，省略时遵循原版默认值。
 * @return 成功解析时返回参数个数，语法错误返回 0。
 */
int CurveFitter::doCustomFit(const QString &equation, const Array &initial)
{
    customFormula = QString();
    customParamCount = getNumParams(equation);
    if (customParamCount == 0)
        return 0;
    auto engine = std::make_shared<QJSEngine>();
    engine->evaluate(
        "var PI=Math.PI; var abs=Math.abs, exp=Math.exp, log=Math.log, sqrt=Math.sqrt, "
        "pow=Math.pow, "
        "sin=Math.sin, cos=Math.cos, tan=Math.tan, asin=Math.asin, acos=Math.acos, atan=Math.atan, "
        "atan2=Math.atan2, floor=Math.floor, ceil=Math.ceil, round=Math.round, minOf=Math.min, "
        "maxOf=Math.max;");
    engine->evaluate(
        "Math.erf=function(x){var A=.05857895078654250866288,B=-.00626245895772819579,"
        "C=-.00299946450696036814,D=.289389696496082416,E=.0539962589851632982,"
        "F=.00508516909930653109,G=.000215969713046142876,H=-.000225663858340491571,I=-3."
        "06833213472529049e-7;"
        "var x2=x*x,k=2/Math.sqrt(Math.PI);if(x2<1e-8)return k*x*(1+x2*(-1/3+x2/10));"
        "var v=x*((k-A)+A*(1+x2*x2*(B+x2*(C+x2*(H+x2*I))))/(1+x2*(D+x2*(E+x2*(F+x2*G)))));"
        "var e=x2>36?1:Math.sqrt(1-Math.exp(-v*v));return x>0?e:-e;};var erf=Math.erf;");
    const QString source = "(function(x,a,b,c,d,e,f){var y=0;" + equation + ";return y;})";
    auto function = engine->evaluate(source);
    if (function.isError() || !function.isCallable())
        return 0;
    customFormula = equation;
    customFunction = [engine, function, count = customParamCount](const Array &params,
                                                                  double x) mutable {
        QJSValueList args;
        args << x;
        for (int i = 0; i < 6; ++i)
            args << (i < count ? params.value(i) : 0.0);
        auto result = function.call(args);
        return result.isError() ? ImageJNumeric::qNaN : result.toNumber();
    };
    macroDefined = true;
    initialParams = initial;
    doFit(CUSTOM);

    return customParamCount;
}
