#include "curvefitterwidgets.h"

#include <atomic>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <QApplication>
#include <QDialog>
#include <QDialogButtonBox>
#include <QEvent>
#include <QFormLayout>
#include <QKeyEvent>
#include <QLabel>
#include <QLineEdit>
#include <QPainter>
#include <QPainterPath>
#include <QPointer>
#include <QSpinBox>
#include <QThread>

namespace {

/**
 * @brief 在 GUI 线程接收 ESC，供优化线程安全读取并清除中止请求。
 */
class EscapeFilter : public QObject
{
public:
    explicit EscapeFilter(QObject *parent)
        : QObject(parent)
    {
    }

    std::atomic<bool> pressed{};

protected:
    bool eventFilter(QObject *, QEvent *event) override
    {
        if (event->type() == QEvent::KeyPress
            && static_cast<QKeyEvent *>(event)->key() == Qt::Key_Escape)
            pressed = true;

        return false;
    }
};

} // namespace

/**
 * @brief 安装 Qt 对话框与线程安全的状态回调，保留原版“取消设置后继续拟合”的行为。
 * @param fitter 要配置的拟合器。
 * @param status 可选状态标签，优化线程通过队列更新它。
 * @param parent 对话框父窗口；设置对话框必须在 GUI 线程调用。
 */
void CurveFitterWidgets::install(CurveFitter &fitter, QLabel *status, QWidget *parent)
{
    if (!qobject_cast<QApplication *>(QCoreApplication::instance()))
        throw std::logic_error("Qt Widgets integration requires QApplication");
    auto filter = std::make_shared<EscapeFilter>(nullptr);
    qApp->installEventFilter(filter.get());
    QPointer<QLabel> label(status);
    fitter.getMinimizer().setStatusCallbacks(
        [label](const QString &text) {
            if (label)
                QMetaObject::invokeMethod(
                    label.data(),
                    [label, text] {
                        if (label)
                            label->setText(text);
                    },
                    Qt::QueuedConnection);
        },
        [filter] { return filter->pressed.exchange(false); });
    const QPointer<QWidget> window(parent);
    fitter.setSettingsCallback([window, &fitter](QVector<double> &params, QVector<double> &,
                                                 int &iterations, int &restarts, double &error) {
        if (QThread::currentThread() != qApp->thread())
            throw std::logic_error("Fit settings must run on the GUI thread");
        QDialog dialog(window);
        dialog.setWindowTitle(
            QCoreApplication::translate("CurveFitter", "Simplex Fitting Options"));
        QFormLayout layout(&dialog);
        layout.addRow(
            new QLabel(QCoreApplication::translate("CurveFitter", "Function name: %1\nFormula: %2")
                           .arg(fitter.getName(), fitter.getFormula()),
                       &dialog));
        QVector<QLineEdit *> fields;
        for (int i = 0; i < params.size(); ++i)
        {
            auto field = new QLineEdit(QString::number(params[i], 'g', 17), &dialog);
            fields << field;
            layout.addRow(QString("Initial %1:").arg(QChar('a' + i)), field);
        }
        QSpinBox maximum, retries;
        maximum.setRange(1, std::numeric_limits<int>::max());
        maximum.setValue(iterations);
        retries.setRange(0, std::numeric_limits<int>::max());
        retries.setValue(restarts);
        QLineEdit exponent(QString::number(-std::log(error) / std::log(10.0), 'g', 17));
        layout.addRow(QCoreApplication::translate("CurveFitter", "Maximum iterations:"), &maximum);
        layout.addRow(QCoreApplication::translate("CurveFitter", "Number of restarts:"), &retries);
        layout.addRow(QCoreApplication::translate("CurveFitter", "Error tolerance [1*10^(-x)]:"),
                      &exponent);
        QDialogButtonBox buttons(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
        layout.addRow(&buttons);
        QObject::connect(&buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
        QObject::connect(&buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
        if (dialog.exec() != QDialog::Accepted)
            return false;
        for (int i = 0; i < params.size(); ++i)
        {
            bool ok{};
            const double value = fields[i]->text().toDouble(&ok);
            if (ok && !std::isnan(value))
                params[i] = value;
        }
        iterations = maximum.value();
        restarts = retries.value();
        bool ok{};
        const double value = exponent.text().toDouble(&ok);
        if (ok)
            error = std::pow(10., -value);

        return true;
    });
}

////////////////////////////////////////////////////////////////////////
// class CurvePlotWidget

CurvePlotWidget::CurvePlotWidget(CurveFitter::PlotData data, QWidget *parent)
    : QWidget(parent)
    , data(std::move(data))
{
    resize(this->data.width, this->data.height);
    setWindowTitle(this->data.title);
}

void CurvePlotWidget::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.fillRect(rect(), Qt::white);
    const QRectF frame(55, 20, width() - 80, height() - 65);
    painter.setPen(Qt::black);
    painter.drawRect(frame);
    painter.drawText(QRectF(0, height() - 35, width(), 25), Qt::AlignCenter, data.xLabel);
    painter.drawText(QPointF(15, 20), data.yLabel);
    const auto position = [&](const QPointF &point) {
        const double dx = data.xMax == data.xMin ? 1 : data.xMax - data.xMin;
        const double dy = data.yMax == data.yMin ? 1 : data.yMax - data.yMin;
        return QPointF(frame.left() + (point.x() - data.xMin) / dx * frame.width(),
                       frame.bottom() - (point.y() - data.yMin) / dy * frame.height());
    };
    painter.save();
    painter.setClipRect(frame);
    QPainterPath curve;
    bool connected{};
    foreach (const QPointF &point, data.curve)
    {
        if (!std::isfinite(point.x()) || !std::isfinite(point.y()))
        {
            connected = false;
            continue;
        }
        if (connected)
            curve.lineTo(position(point));
        else
            curve.moveTo(position(point));
        connected = true;
    }
    painter.setPen(QColor(0x0000FF));
    painter.drawPath(curve);
    painter.setPen(QColor(0xFF0000));
    painter.setBrush(QColor(0xFF0000));
    foreach (const QPointF &point, data.data)
        if (std::isfinite(point.x()) && std::isfinite(point.y()))
            painter.drawEllipse(position(point), 2.5, 2.5);
    painter.restore();
    painter.setPen(QColor(0x0000FF));
    painter.drawText(frame.adjusted(10, 10, -10, -10), Qt::AlignLeft | Qt::AlignTop, data.legend);
}
