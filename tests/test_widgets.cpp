#include "curvefitterwidgets.h"

#include <QApplication>
#include <QDialog>
#include <QImage>
#include <QLabel>
#include <QTest>
#include <QTimer>

/**
 * @brief 验证原版设置取消行为及 Qt 曲线绘图适配。
 */
class WidgetsTest : public QObject
{
    Q_OBJECT

private slots:
    void settingsAndPlot();
};

void WidgetsTest::settingsAndPlot()
{
    CurveFitter fitter({1, 2, 3}, {2, 4, 8});
    QLabel status;
    CurveFitterWidgets::install(fitter, &status);
    bool displayed{};
    QTimer::singleShot(0, [&] {
        auto dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget());
        if (dialog)
        {
            displayed = true;
            dialog->reject();
        }
    });
    fitter.doFit(CurveFitter::EXPONENTIAL, true);
    QVERIFY(displayed);
    QCOMPARE(fitter.getStatus(), Minimizer::SUCCESS);
    CurvePlotWidget plot(fitter.getPlot());
    const QImage image = plot.grab().toImage();
    QVERIFY(!image.isNull());
    int red{}, blue{};
    for (int y = 0; y < image.height(); ++y)
        for (int x = 0; x < image.width(); ++x)
        {
            const auto color = image.pixelColor(x, y);
            if (color.red() > 200 && color.blue() < 50)
                ++red;
            if (color.blue() > 200 && color.red() < 50)
                ++blue;
        }
    QVERIFY(red > 0);
    QVERIFY(blue > 0);
}

QTEST_MAIN(WidgetsTest)
#include "test_widgets.moc"
