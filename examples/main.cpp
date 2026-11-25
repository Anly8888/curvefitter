#include "curvefitter.h"
#include <QTextStream>

int main()
{
    CurveFitter fitter({0, 1, 2, 3, 4}, {2, 5, 8, 11, 14});
    fitter.doFit(CurveFitter::STRAIGHT_LINE);
    QTextStream(stdout) << fitter.getResultString() << '\n';
    return fitter.getStatus() == CurveFitter::SUCCESS ? 0 : 1;
}
