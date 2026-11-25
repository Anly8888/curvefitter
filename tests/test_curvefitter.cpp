#include "curvefitter.h"

#include <cmath>
#include <thread>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QTest>

/**
 * @brief 对照未修改的 ImageJ Java 实现，验证 Qt 移植的模型和优化器行为。
 */
class CurveFitterTest : public QObject
{
    Q_OBJECT

private slots:
    void javaReference_data();
    void javaReference();
    void customFunction();
    void formula();
    void parallelSearch();
    void resolutionAndExtraElements();
    void plot();
    void settings();
};

namespace {

QVector<double> numbers(const QJsonValue &value)
{
    QVector<double> result;
    foreach (const QJsonValue &item, value.toArray())
        result << (item.isNull() ? std::numeric_limits<double>::quiet_NaN() : item.toDouble());

    return result;
}

void compareNumber(double actual, double expected, double tolerance = 1e-5)
{
    if (std::isnan(expected))
        QVERIFY(std::isnan(actual));
    else
        QVERIFY2(std::abs(actual - expected) <= tolerance * (1 + std::abs(expected)),
                 qPrintable(QString("actual=%1 expected=%2")
                                .arg(actual, 0, 'g', 17)
                                .arg(expected, 0, 'g', 17)));
}

} // namespace

void CurveFitterTest::javaReference_data()
{
    QTest::addColumn<QByteArray>("record");
    QFile file(QStringLiteral(REFERENCE_PATH));
    QVERIFY2(file.open(QIODevice::ReadOnly), qPrintable(file.errorString()));
    while (!file.atEnd())
    {
        const QByteArray line = file.readLine().trimmed();
        const auto object = QJsonDocument::fromJson(line).object();
        QTest::newRow(qPrintable(object["name"].toString())) << line;
    }
}

void CurveFitterTest::javaReference()
{
    QFETCH(QByteArray, record);
    const auto reference = QJsonDocument::fromJson(record).object();
    QVector<double> actual;
    if (reference.contains("minimizer"))
    {
        const int mode = reference["minimizer"].toInt();
        Minimizer minimizer;
        minimizer.setMaximumThreads(1);
        minimizer.setRandomSeed(31);
        int calls{};
        minimizer.setFunction(
            [&](QVector<double> &p, double) {
                if (mode == 4 && ++calls > 100)
                    return std::numeric_limits<double>::quiet_NaN();
                if (mode == 1)
                    return std::numeric_limits<double>::quiet_NaN();
                if (mode == 2)
                {
                    minimizer.abort();
                    return p[0] * p[0];
                }

                return (p[0] - 2) * (p[0] - 2) + (p[1] + 3) * (p[1] + 3)
                       + (mode == 5 ? minimizer.getCompletedMinimizations() : 0);
            },
            2);
        if (mode == 3)
            minimizer.setMaxIterations(1);
        if (mode == 5)
            minimizer.setMaxRestarts(1);
        minimizer.setMaxError(1e-10, 1e-12);
        QCOMPARE(mode == 4 ? minimizer.minimizeOnce({0, 0}, {1, 1})
                           : minimizer.minimize({0, 0}, {1, 1}),
                 reference["status"].toInt());
        QCOMPARE(minimizer.getIterations(), reference["iterations"].toInt());
        QCOMPARE(minimizer.getCompletedMinimizations(), reference["completed"].toInt());
        actual = minimizer.getParams();
    }
    else
    {
        CurveFitter fitter(numbers(reference["x"]), numbers(reference["y"]));
        fitter.getMinimizer().setMaximumThreads(1);
        fitter.getMinimizer().setRandomSeed(0);
        if (!reference["initial"].isNull())
            fitter.setInitialParameters(numbers(reference["initial"]));
        if (!reference["weights"].isNull())
            fitter.setWeights(numbers(reference["weights"]));
        fitter.doFit(reference["type"].toInt());
        QCOMPARE(fitter.getStatus(), reference["status"].toInt());
        actual = fitter.getParams();
        compareNumber(fitter.getSumResidualsSqr(), reference["sse"].toDouble());
        compareNumber(fitter.getRSquared(), reference["r2"].toDouble());
        compareNumber(fitter.getSD(), reference["sd"].toDouble());
        compareNumber(fitter.getFitGoodness(), reference["goodness"].toDouble());
        const auto expected = numbers(reference["residuals"]);
        const auto residuals = fitter.getResiduals();
        QCOMPARE(residuals.size(), expected.size());
        for (int i = 0; i < expected.size(); ++i)
            compareNumber(residuals[i], expected[i]);
    }
    const auto expected = numbers(reference["params"]);
    QCOMPARE(actual.size(), expected.size());
    for (int i = 0; i < expected.size(); ++i)
        compareNumber(actual[i], expected[i]);
}

void CurveFitterTest::customFunction()
{
    CurveFitter fitter({0, 1, 2, 3, 4}, {1, 4, 9, 16, 25});
    fitter.setOffsetMultiplySlopeParams(0, -1, 1);
    fitter.doCustomFit([](const QVector<double> &p,
                          double x) { return p[0] + p[1] * x + p[2] * x * x; },
                       3, "a+b*x+c*x*x", {0, 0, 1});
    QCOMPARE(fitter.getStatus(), Minimizer::SUCCESS);
    compareNumber(fitter.getParams()[0], 1);
    compareNumber(fitter.getParams()[1], 2);
    compareNumber(fitter.getParams()[2], 1);
}

void CurveFitterTest::formula()
{
    CurveFitter fitter({0, 1, 2, 3}, {2, 5, 8, 11});
    QCOMPARE(fitter.doCustomFit("y=a+b*x", {1, 1}), 2);
    QCOMPARE(fitter.getStatus(), Minimizer::SUCCESS);
    compareNumber(fitter.f(4), 14, 1e-4);
    QCOMPARE(fitter.doCustomFit("y = a +", {1}), 0);
    QCOMPARE(CurveFitter::getNumParams("y=a+Math.exp(b*x)"), 2);
}

void CurveFitterTest::parallelSearch()
{
    Minimizer minimizer;
    minimizer.setMaximumThreads(2);
    minimizer.setFunction([](QVector<double> &p,
                             double) { return std::pow(p[0] - 2, 2) + std::pow(p[1] + 3, 2); },
                          2);
    minimizer.setMaxError(1e-10, 1e-12);
    QCOMPARE(minimizer.minimize({0, 0}, {1, 1}), Minimizer::SUCCESS);
    QVERIFY(minimizer.getCompletedMinimizations() >= 2);
    compareNumber(minimizer.getParams()[0], 2);
    compareNumber(minimizer.getParams()[1], -3);
}

void CurveFitterTest::resolutionAndExtraElements()
{
    Minimizer minimizer;
    minimizer.setMaximumThreads(1);
    minimizer.setExtraArrayElements(2);
    minimizer.setParamResolutions({1e-5});
    minimizer.setFunction(
        [](QVector<double> &p, double) {
            p[2] = 42;
            p[3] = 81;
            return (p[0] - 2) * (p[0] - 2);
        },
        1);
    QCOMPARE(minimizer.minimize({0}, {1}), Minimizer::SUCCESS);
    QCOMPARE(minimizer.getParams().size(), 4);
    QCOMPARE(minimizer.getParams()[2], 42.0);
    QCOMPARE(minimizer.getParams()[3], 81.0);
    compareNumber(minimizer.getFunctionValue(), 0);
}

void CurveFitterTest::plot()
{
    CurveFitter fitter({0, 1, 2}, {2, 5, 8});
    fitter.doFit(CurveFitter::STRAIGHT_LINE);
    const auto data = fitter.getPlot(256);
    QCOMPARE(data.curve.size(), 256);
    QCOMPARE(data.curve.first(), QPointF(0, 2));
    QCOMPARE(data.curve.last(), QPointF(255, 767));
    QCOMPARE(data.data.size(), 3);
    QVERIFY(fitter.getResultString().contains("Number of iterations: 1"));
}

void CurveFitterTest::settings()
{
    CurveFitter fitter({1, 2, 3}, {2, 4, 8});
    bool called{};
    fitter.setSettingsCallback([&](QVector<double> &, QVector<double> &, int &, int &, double &) {
        called = true;
        return false;
    });
    fitter.doFit(CurveFitter::EXPONENTIAL, true);
    QVERIFY(called);
    QCOMPARE(fitter.getStatus(), Minimizer::SUCCESS);
}

QTEST_GUILESS_MAIN(CurveFitterTest)
#include "test_curvefitter.moc"
