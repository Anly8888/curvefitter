#include "curvefitter.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
static int checks = 0;

void check(bool ok, const char* message)
{
    ++checks;
    if (!ok)
        throw std::runtime_error(message);
}

void near(double a, double b, double tol = 1e-6)
{
    check(std::isfinite(a) && std::abs(a - b) < tol, "Numerical mismatch");
}

int main()
{
    try
    {
        for (int type = 0; type < 25; ++type)
        {
            QVector<double> p;
            switch (type)
            {
            case 0:
                p = { 2, 3 };
                break;
            case 1:
                p = { 2, 3, .2 };
                break;
            case 2:
                p = { 2, 3, .2, -.1 };
                break;
            case 3:
                p = { 2, 3, .2, -.1, .05 };
                break;
            case 4:
            case 15:
                p = { 2, .3 };
                break;
            case 5:
            case 16:
                p = { 2, 1.3 };
                break;
            case 6:
                p = { 2, 3 };
                break;
            case 7:
                p = { 8, 2, 3, 1 };
                break;
            case 8:
                p = { .1, 2, 2, 1 };
                break;
            case 9:
                p = { 2, 3, -1 };
                break;
            case 10:
            case 14:
                p = { 0, 2, 3, 8 };
                break;
            case 11:
            case 13:
                p = { 2, .7, 1 };
                break;
            case 12:
                p = { 1, 4, 2, .7 };
                break;
            case 21:
                p = { 3, 2, .7 };
                break;
            case 22:
                p = { 3, .7 };
                break;
            case 23:
                p = { 3, .7, 1.5 };
                break;
            case 24:
                p = { 2, 3, 2, .8 };
                break;
            default:
                p.fill(.01, CurveFitter::getNumParams(type));
                p[0] = 2;
                p[1] = 1;
            }
            QVector<double> x, y;
            for (int i = 0; i < 60; ++i)
            {
                double xx = .2 + 4. * i / 59;
                x << xx;
                y << CurveFitter::f(type, p, xx);
            }
            CurveFitter cf(x, y);
            cf.doFit(type);
            if (cf.getStatus() != CurveFitter::SUCCESS || cf.getSumResidualsSqr() > 1e-5)
            {
                std::cerr << "Default fit failed type " << type << " "
                          << cf.getStatusString().toStdString() << " SSE "
                          << cf.getSumResidualsSqr() << '\n';
                return 1;
            }
            ++checks;
            cf.setInitialParameters(p);
            cf.doFit(type);
            check(cf.getStatus() == 0, "Seeded fit failed");
            for (int i = 0; i < x.size(); ++i)
                near(cf.f(cf.getParams(), x[i]), y[i], 1e-4);
        }
        CurveFitter line({ 0, 1, 2 }, { 1, 3, 100 });
        line.setWeights({ 1, 1, 0 });
        line.doFit(0);
        near(line.getParams()[0], 1);
        near(line.getParams()[1], 2);
        CurveFitter clean({ 0, NAN, 1 }, { 1, 99, 3 });
        clean.setWeights({ 1, 0, 1 });
        clean.doFit(0);
        near(clean.getParams()[1], 2);
        CurveFitter custom({ 0, 1, 2, 3 }, { 1, 4, 9, 16 });
        custom.doCustomFit([](const QVector<double>& p, double x)
            { return p[0] + p[1] * x + p[2] * x * x; }, 3, "a+b*x+c*x*x", { .5, 1, 1 });
        check(custom.getStatus() == 0, "Custom fit failed");
        near(custom.getParams()[1], 2, 1e-4);
        CurveFitter invalid({ 0, 0 }, { 1, 2 });
        invalid.doFit(0);
        check(invalid.getStatus() == 1, "Singular input accepted");
        CurveFitter negative({ 1, 2, 3 }, { -2, -4, -8 });
        negative.doFit(15);
        near(negative.getParams()[0], -1);
        near(negative.getParams()[1], std::log(2.));
        CurveFitter origin({ 0, 1, 2, 3 }, { 0, 2, 8, 18 });
        origin.doFit(16);
        near(origin.getParams()[1], 2);
        CurveFitter domain({ -1, 1, 2 }, { 2, 3, 4 });
        domain.doFit(6);
        check(domain.getStatus() == 1, "Log domain accepted");
        line.setMaxIterations(1);
        line.doCustomFit([](const QVector<double>& p, double x)
            { return p[0] * std::exp(p[1] * x); }, 2, "a*exp(b*x)", { 1, .1 });
        check(line.getStatus() == 2, "Iteration limit ignored");
        check(CurveFitter::fitList().size() == 25, "Missing model metadata");

        // Independent reference from SciPy least_squares (not the built-in evaluator).
        CurveFitter noisy({ 0.0, 0.3, 0.6, 0.8999999999999999, 1.2, 1.5, 1.7999999999999998, 2.1,
                              2.4, 2.6999999999999997 },
            { 2.1, 2.1349937031587514, 2.5924983006428093, 2.89665882912068, 3.162148804385787,
                3.7242376007810183, 4.068866421287775, 4.722733953562184, 5.123392946846236,
                5.909359102131046 });
        noisy.doFit(CurveFitter::EXPONENTIAL);
        near(noisy.getParams()[0], 2.0108545150872716, 1e-5);
        near(noisy.getParams()[1], 0.3973396954541083, 1e-5);
        near(CurveFitter::f(CurveFitter::GAUSSIAN, { 1, 5, 2, 3 }, 2), 5);
        near(CurveFitter::f(CurveFitter::RODBARD, { 8, 2, 3, 2 }, 3), 5);
        near(CurveFitter::f(CurveFitter::GAMMA_VARIATE, { 1, 2, 2, 3 }, 0), 0);
        std::cout << checks << " checks passed\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
