#pragma once

#include "curvefitter.h"

#include <QWidget>

class QLabel;

/**
 * @brief 将 ImageJ 设置对话框、状态行及 ESC 中止映射到 Qt Widgets。
 */
namespace CurveFitterWidgets {

void install(CurveFitter &fitter, QLabel *status = nullptr, QWidget *parent = nullptr);

}

/**
 * @brief 绘制 CurveFitter 的拟合曲线、原始散点和参数图例。
 */
class CurvePlotWidget : public QWidget
{
public:
    explicit CurvePlotWidget(CurveFitter::PlotData data, QWidget *parent = nullptr);

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    CurveFitter::PlotData data;
};
