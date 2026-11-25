#include "curvefitter.h"
#include "numericcompat.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <QElapsedTimer>
#include <QRegularExpression>

using namespace ImageJNumeric;

////////////////////////////////////////////////////////////////////////
// class CurveFitter

QStringList CurveFitter::fitList()
{
    return {"Straight Line",
            "2nd Degree Polynomial",
            "3rd Degree Polynomial",
            "4th Degree Polynomial",
            "Exponential",
            "Power",
            "Log",
            "Rodbard",
            "Gamma Variate",
            "y = a+b*ln(x-c)",
            "Rodbard (NIH Image)",
            "Exponential with Offset",
            "Gaussian",
            "Exponential Recovery",
            "Inverse Rodbard",
            "Exponential (linear regression)",
            "Power (linear regression)",
            "5th Degree Polynomial",
            "6th Degree Polynomial",
            "7th Degree Polynomial",
            "8th Degree Polynomial",
            "Gaussian (no offset)",
            "Exponential Recovery (no offset)",
            "Chapman-Richards",
            "Error Function"};
}

QStringList CurveFitter::formulas()
{
    return {"y = a+bx",
            "y = a+bx+cx^2",
            "y = a+bx+cx^2+dx^3",
            "y = a+bx+cx^2+dx^3+ex^4",
            "y = a*exp(bx)",
            "y = a*x^b",
            "y = a*ln(bx)",
            "y = d+(a-d)/(1+(x/c)^b)",
            "y = b*(x-a)^c*exp(-(x-a)/d)",
            "y = a+b*ln(x-c)",
            "x = d+(a-d)/(1+(y/c)^b) [y = c*((x-a)/(d-x))^(1/b)]",
            "y = a*exp(-bx) + c",
            "y = a + (b-a)*exp(-(x-c)*(x-c)/(2*d*d))",
            "y = a*(1-exp(-b*x)) + c",
            "y = c*((x-a)/(d-x))^(1/b)",
            "y = a*exp(bx)",
            "y = a*x^b",
            "y = a+bx+cx^2+dx^3+ex^4+fx^5",
            "y = a+bx+cx^2+dx^3+ex^4+fx^5+gx^6",
            "y = a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7",
            "y = a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7+ix^8",
            "y = a*exp(-(x-b)*(x-b)/(2*c*c))",
            "y = a*(1-exp(-b*x))",
            "y = a*(1-exp(-b*x))^c",
            "y = a+b*erf((x-c)/d)"};
}

QStringList CurveFitter::macroCodes()
{
    return {"y = a+x*b",
            "y = a+x*(b+x*c)",
            "y = a+x*(b+x*(c+x*d))",
            "y = a+x*(b+x*(c+x*(d+x*e)))",
            "y = a*Math.exp(b*x)",
            "y = a*Math.pow(x,b)",
            "y = a*Math.log(b*x)",
            "y = d+(a-d)/(1+Math.pow(x/c,b))",
            "y = b*Math.pow(x-a,c)*Math.exp(-(x-a)/d)",
            "y = a+b*Math.log(x-c)",
            "y = c*Math.pow((x-a)/(d-x),1/b)",
            "y = a*Math.exp(-b*x)+c",
            "y = a+(b-a)*Math.exp(-(x-c)*(x-c)/(2*d*d))",
            "y = a*(1-Math.exp(-b*x))+c",
            "y = c*Math.pow((x-a)/(d-x),1/b)",
            "y = a*Math.exp(b*x)",
            "y = a*Math.pow(x,b)",
            "y = a+x*(b+x*(c+x*(d+x*(e+x*f))))",
            "y = a+x*(b+x*(c+x*(d+x*(e+x*(f+x*g)))))",
            "y = a+x*(b+x*(c+x*(d+x*(e+x*(f+x*(g+x*h))))))",
            "y = a+x*(b+x*(c+x*(d+x*(e+x*(f+x*(g+x*(h+x*i)))))))",
            "y = a*Math.exp(-(x-b)*(x-b)/(2*c*c))",
            "y = a*(1-Math.exp(-b*x))",
            "y = a*Math.pow(1-Math.exp(-b*x),c)",
            "y = a+b*Math.erf((x-c)/d)"};
}

QVector<int> CurveFitter::sortedTypes()
{
    return {0,  1,  2, 3, 17, 18, 19, 20, 5,  16, 4, 15, 11,
            13, 22, 6, 9, 12, 21, 24, 7,  10, 14, 8, 23};
}

QStringList CurveFitter::getSortedFitList()
{
    QStringList result;
    foreach (int type, sortedTypes())
        result << fitList()[type];
    return result;
}

int CurveFitter::getFitCode(const QString &name)
{
    return fitList().indexOf(name);
}

CurveFitter::CurveFitter(const Array &x, const Array &y)
{
    if (x.size() != y.size())
        throw std::invalid_argument("x/y size mismatch");
    for (int i = 0; i < x.size(); ++i)
        if (!std::isnan(x[i] + y[i]))
        {
            xData << x[i];
            yData << y[i];
        }
    numPoints = xData.size();
}

CurveFitter::Array CurveFitter::getParams()
{
    return finalParams.isEmpty() ? minimizer.getParams() : finalParams;
}

QString CurveFitter::getStatusString() const
{
    return errorString.isNull() ? Minimizer::statusStrings().value(getStatus()) : errorString;
}

QString CurveFitter::getName()
{
    if (fitType == CUSTOM)
        return tr("User-defined");
    if (fitType == GAUSSIAN_INTERNAL)
        fitType = GAUSSIAN;
    else if (fitType == RODBARD_INTERNAL)
        fitType = RODBARD;

    return fitList().value(fitType);
}

QString CurveFitter::getFormula()
{
    if (fitType == CUSTOM)
        return customFormula;
    if (fitType == GAUSSIAN_INTERNAL)
        fitType = GAUSSIAN;
    else if (fitType == RODBARD_INTERNAL)
        fitType = RODBARD;

    return formulas().value(fitType);
}

QString CurveFitter::getMacroCode()
{
    if (fitType == CUSTOM)
        return customFormula;
    if (fitType == GAUSSIAN_INTERNAL)
        fitType = GAUSSIAN;
    else if (fitType == RODBARD_INTERNAL)
        fitType = RODBARD;

    return macroCodes().value(fitType);
}

double CurveFitter::f(double x)
{
    if (finalParams.isEmpty())
        finalParams = minimizer.getParams();
    return f(finalParams, x);
}

double CurveFitter::f(const Array &p, double x)
{
    return fitType == CUSTOM ? customFunction(p, x) : f(fitType, p, x);
}

void CurveFitter::doCustomFit(UserFunction function,
                              int count,
                              const QString &formula,
                              const Array &initial,
                              const Array &variations,
                              bool showSettings)
{
    customFunction = std::move(function);
    customParamCount = count;
    initialParams = initial;
    initialParamVariations = variations;
    customFormula = formula.isNull() ? tr("(defined in plugin)") : formula;
    doFit(CUSTOM, showSettings);
}

void CurveFitter::settingsDialog()
{
    if (!settingsCallback)
        throw std::logic_error("Install a Qt settings callback before showSettings=true");
    if (initialParamVariations.isEmpty())
        initialParamVariations = Array(numParams);
    int iterations = minimizer.getMaxIterations(), restarts = minimizer.getMaxRestarts();
    double error = maxRelError;
    const Array previous = initialParams;
    if (!settingsCallback(initialParams, initialParamVariations, iterations, restarts, error))
    {
        initialParams = previous;
        return;
    }
    for (int i = 0; i < numParams; ++i)
        if (!std::isnan(initialParams[i]))
            initialParamVariations[i] = jmax(.01 * initialParams[i],
                                             .001 * initialParamVariations[i]);
    if (iterations > 0)
        minimizer.setMaxIterations(iterations);
    if (restarts >= 0)
        minimizer.setMaxRestarts(restarts);
    setMaxError(error);
}

void CurveFitter::doFit(int fitType, bool showSettings)
{
    if (!(fitType >= STRAIGHT_LINE && fitType < fitList().size() || fitType == CUSTOM))
        throw std::invalid_argument("Invalid fit type");
    if (fitType == CUSTOM && !customFunction)
        throw std::invalid_argument("No custom formula!");
    this->fitType = fitType;
    if (isModifiedFitType(
            fitType)) // these fits don't use the original data and a different fit type (this->fitType)
        if (!prepareModifiedFitType(fitType))
            return;
    numParams = getNumParams();
    if (fitType != CUSTOM)
        getOffsetAndFactorParams();

    calculateSumYandY2(); // sumY, sumY2 needed for regression, abs Error; R, goodness of modified
                          // fit functions
    QElapsedTimer elapsed;
    elapsed.start();
    if (this->fitType == STRAIGHT_LINE)
    {                                 // no minimizer needed
        finalParams = Array{0, 0, 0}; // (does not work with 0 initial slope)
        doRegression(finalParams); // fit by regression; save sum(Residuals^2) as last array element
        linearRegressionUsed = true;
    }
    else
    { // use simplex minimizer
        minimizer.setFunction([this](Array &params, double x) { return userFunction(params, x); },
                              numParams - numRegressionParams);
        minimizer.setExtraArrayElements(
            numRegressionParams); // reserve space for extra parameters if minimizer has fewer paramters
        if (macroDefined)
            minimizer.setMaximumThreads(1); // macro interpreter does not allow multithreading
        if (!makeInitialParamsAndVariations(fitType)) // also includes some data checking
            return;                                   // initialization failure
        if (showSettings)
            settingsDialog();
        if (numRegressionParams > 0)
            modifyInitialParamsAndVariations();
        else
        {
            minimizerInitialParams = initialParams;
            minimizerInitialParamVariations = initialParamVariations;
        }
        elapsed.restart();
        // The maximum absolute error of the fit must be specified in case the
        // fit function fits perfectly, i.e. the sume of residuals approaches 0.
        // In such a case, the maximum relative error is meaningless and the
        // minimizer would run until it reaches the maximum iteration count.
        double maxAbsError = jmin(1e-6, maxRelError) * std::sqrt(sumY2);
        minimizer.setMaxError(maxRelError, maxAbsError);

        // m i n i m i z a t i o n	of squared residuals
        minimizerStatus = minimizer.minimize(minimizerInitialParams,
                                             minimizerInitialParamVariations);
        finalParams = minimizer.getParams();
        if (numRegressionParams > 0)
            minimizerParamsToFullParams(finalParams, false);
    }
    if (isModifiedFitType(fitType)) // params of actual fit to user params
        postProcessModifiedFitType(fitType);
    if (fitType == ERF)             // make it nicer
        if (finalParams[3] < 0)
        {                           // y = a+b*erf((x-c)/d) with negative d: invert b instead
            finalParams[1] = -finalParams[1];
            finalParams[3] = -finalParams[3];
        }

    switch (fitType)
    {              // postprocessing for nicer output:
    case GAUSSIAN: // Gaussians: width (std deviation) should be >0
        finalParams[3] = std::abs(finalParams[3]);
        break;
    case GAUSSIAN_NOOFFSET:
        finalParams[2] = std::abs(finalParams[2]);
        break;
    }
    time = elapsed.elapsed();
}

/**
 * @brief For improved fitting performance when using a custom fit formula, one may
 *	specify parameters that can be calculated directly by linear regression.
 *	For values not used, set the index to -1
 *
 * @param offsetParam  Index of a parameter that is a pure offset:
 *					   E.g. '0' if	f(p0, p1, p2...) = p0 + function(p1, p2, ...).
 * @param multiplyParam	 Index of a parameter that is purely multiplicative.
 *					   E.g. multiplyParams=1 if f(p0, p1, p2, p3...) can be expressed as
 *					   p1*func(p0, p2, p3, ...) or p0 +p1*func(p0, p2, p3, ...) with '0' being
 *					   the offsetparam.
 * @param slopeParam   Index of a parameter that is multiplied with x and then summed to the
 * function. E.g. '1' for f(p0, p1, p2, p3...) = p1*x + func(p0, p2, p3, ...) Only one,
 * multiplyParam and slopeParam can be used (ie.e, the other should be set to -1)
 */
void CurveFitter::setOffsetMultiplySlopeParams(int offsetParam, int multiplyParam, int slopeParam)
{
    if (multiplyParam >= 0 && slopeParam >= 0)
        throw std::invalid_argument(
            "CurveFitter: only one of multiplyParam and slopeParam may be given (i.e., >=0)");
    this->offsetParam = offsetParam;
    hasSlopeParam = slopeParam >= 0;
    factorParam = hasSlopeParam ? slopeParam : multiplyParam;
    numRegressionParams = 0;
    if (factorParam >= 0 && factorParam == offsetParam)
        throw std::invalid_argument(
            "CurveFitter: offsetParam and slopeParam/factorParam must be different");
    if (offsetParam >= 0)
        numRegressionParams++;
    if (factorParam >= 0)
        numRegressionParams++;
}

/**
 * @brief Get number of parameters for current fit formula
 *	Do not use before 'doFit', because the fit function would be undefined.	 */
int CurveFitter::getNumParams()
{
    if (fitType == CUSTOM)
        return customParamCount;
    else
        return getNumParams(fitType);
}

int CurveFitter::getNumParams(int fitType)
{
    switch (fitType)
    {
    case STRAIGHT_LINE:
        return 2;
    case POLY2:
        return 3;
    case POLY3:
        return 4;
    case POLY4:
        return 5;
    case POLY5:
        return 6;
    case POLY6:
        return 7;
    case POLY7:
        return 8;
    case POLY8:
        return 9;
    case EXPONENTIAL:
    case EXP_REGRESSION:
        return 2;
    case POWER:
    case POWER_REGRESSION:
        return 2;
    case EXP_RECOVERY_NOOFFSET:
        return 2;
    case LOG:
        return 2;
    case LOG2:
        return 3;
    case GAUSSIAN_NOOFFSET:
        return 3;
    case EXP_RECOVERY:
        return 3;
    case CHAPMAN:
        return 3;
    case EXP_WITH_OFFSET:
        return 3;
    case RODBARD:
    case RODBARD2:
    case INV_RODBARD:
    case RODBARD_INTERNAL:
        return 4;
    case GAMMA_VARIATE:
        return 4;
    case GAUSSIAN:
    case GAUSSIAN_INTERNAL:
        return 4;
    case ERF:
        return 4;
    }
    return 0;
}

double CurveFitter::f(int fitType, const Array &p, double x)
{
    switch (fitType)
    {
    case STRAIGHT_LINE:
        return p[0] + x * p[1];
    case POLY2:
        return p[0] + x * (p[1] + x * p[2]);
    case POLY3:
        return p[0] + x * (p[1] + x * (p[2] + x * p[3]));
    case POLY4:
        return p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * p[4])));
    case POLY5:
        return p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * p[5]))));
    case POLY6:
        return p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * (p[5] + x * p[6])))));
    case POLY7:
        return p[0]
               + x
                     * (p[1]
                        + x * (p[2] + x * (p[3] + x * (p[4] + x * (p[5] + x * (p[6] + x * p[7]))))));
    case POLY8:
        return p[0]
               + x
                     * (p[1]
                        + x
                              * (p[2]
                                 + x
                                       * (p[3]
                                          + x
                                                * (p[4]
                                                   + x
                                                         * (p[5]
                                                            + x * (p[6] + x * (p[7] + x * p[8])))))));
    case EXPONENTIAL:
    case EXP_REGRESSION:
        return p[0] * std::exp(p[1] * x);
    case EXP_WITH_OFFSET:
        return p[0] * std::exp(-p[1] * x) + p[2];
    case EXP_RECOVERY:
        return p[0] * (1 - std::exp(-p[1] * x)) + p[2];
    case EXP_RECOVERY_NOOFFSET:
        return p[0] * (1 - std::exp(-p[1] * x));
    case CHAPMAN: // a*(1-exp(-b*x))^c
    {
        double value = p[0] * (std::pow((1 - (std::exp(-p[1] * x))), p[2]));
        //	Log.e("test", "values = " + value);
        return value;
    }
    case GAUSSIAN:
        return p[0] + (p[1] - p[0]) * std::exp(-(x - p[2]) * (x - p[2]) / (2.0 * p[3] * p[3]));
    case GAUSSIAN_INTERNAL: // replaces GAUSSIAN for the fitting process
        return p[0] + p[1] * std::exp(-(x - p[2]) * (x - p[2]) / (2.0 * p[3] * p[3]));
    case GAUSSIAN_NOOFFSET:
        return p[0] * std::exp(-(x - p[1]) * (x - p[1]) / (2.0 * p[2] * p[2]));
    case POWER: // ax^b
    case POWER_REGRESSION:
        return p[0] * std::pow(x, p[1]);
    case LOG:
        if (x == 0.0)
            return -1000 * p[0];
        return p[0] * std::log(p[1] * x);
    case RODBARD:
    {                                         // d+(a-d)/(1+(x/c)^b)
        double ex = std::pow(x / p[2], p[1]); //(x/c)^b
        return p[3] + (p[0] - p[3]) / (1.0 + ex);
    }
    case RODBARD_INTERNAL:
    { // d+a/(1+(x/c)^b) , replaces RODBARD of the fitting process
        double ex = std::pow(x / p[2], p[1]); //(x/c)^b
        return p[3] + p[0] / (1.0 + ex);
    }
    case GAMMA_VARIATE:
    { // b*(x-a)^c*exp(-(x-a)/d)
        if (p[0] >= x)
            return 0.0;
        if (p[1] <= 0)
            return qNaN;
        if (p[2] <= 0)
            return qNaN;
        if (p[3] <= 0)
            return qNaN;

        double pw = std::pow((x - p[0]), p[2]);
        double e = std::exp((-(x - p[0])) / p[3]);
        return p[1] * pw * e;
    }
    case LOG2:
    {
        double tmp = x - p[2];
        if (tmp <= 0)
            return qNaN;
        return p[0] + p[1] * std::log(tmp);
    }
    case INV_RODBARD:    // c*((x-a)/(d-x))^(1/b), the inverse Rodbard function
    case RODBARD2:
    {                    // also used after the 'Rodbard NIH Image' fit
        double y;
        if (p[3] - x < 2 * std::numeric_limits<double>::denorm_min()
            || x < p[0]) // avoid x>=d (singularity) and x<a (negative exponent)
            y = fitType == INV_RODBARD ? qNaN : 0.0;
        else
        {
            y = (x - p[0]) / (p[3] - x); //(x-a)/(d-x) = ( (a-d)/(x-d) -1 )
            y = std::pow(y, 1.0 / p[1]); // y=y**(1/b)
            y = y * p[2];
        }
        return y;
    }
    case ERF:
        return p[0] + p[1] * imagejErf((x - p[2]) / p[3]); // y=a+b*erf((x-c)/d)
    default:
        return 0.0;
    }
}

/**
 * @brief Returns residuals array, i.e., differences between data and curve.
 *	The residuals are with respect to the real data, also for fit types where the data are
 *	modified before fitting (power&exp fit by linear regression, 'Rodbard NIH Image' ).
 *	This is in contrast to sum of squared residuals, which is for the fit that was actually done.
 */
CurveFitter::Array CurveFitter::getResiduals()
{
    Array params = getParams();
    Array residuals = Array(xData.size());
    for (int i = 0; i < xData.size(); i++)
        residuals[i] = yData.constData()[i] - f(params, xData.constData()[i]);
    return residuals;
}

/**
 * @brief Returns the sum of the residuals (may be NaN if the minimizer could not start properly
 *	i.e., if getStatus() returns Minimizer.INITILIZATION_FAILURE).
 *  If weights have been specified, each of the residuals is multiplied by the corresponding
 *  weight before summing.
 */
double CurveFitter::getSumResidualsSqr()
{
    return getParams()[numParams]; // value is stored as last element by the minimizer
}

/**
 * @brief Returns the standard deviation of the residuals.
 *	Here, the standard deviation is defined here as the root-mean-square of the residuals
 *	times sqrt(n/(n-1)); where n is the number of points.
 *  If weights are provided, the standard deviation does not take the weights into account.
 *  With weights, the standard deviation and getSumResidualsSqr (which uses weights)
 *  are not related the usual way.
 */
double CurveFitter::getSD()
{
    Array residuals = getResiduals();
    int n = residuals.size();
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += residuals[i];
        sum2 += residuals[i] * residuals[i];
    }
    double stdDev = (sum2 - sum * sum / n); // sum of squared residuals
    return std::sqrt(stdDev / (n - 1.0));
}

/**
 * @brief Returns R^2, where 1.0 is best. For unweighted data,
    <pre>
     r^2 = 1 - SSE/SSD

     where:  SSE = sum of the squared errors
             SSD = sum of the squared deviations about the mean.
    </pre>
     *	For power, exp by linear regression and 'Rodbard NIH Image', this is calculated for the
     *	fit actually done, not for the residuals of the original data.
     */
double CurveFitter::getRSquared()
{
    if (std::isnan(sumY))
        calculateSumYandY2();
    double sumMeanDiffSqr = sumY2 - sumY * sumY / sumWeights;
    double rSquared = 0.0;
    if (sumMeanDiffSqr > 0.0)
        rSquared = 1.0 - getSumResidualsSqr() / sumMeanDiffSqr;
    return rSquared;
}

/**
 * @brief Get a measure of "goodness of fit" where 1.0 is best.
 *	Approaches R^2 if the number of points is much larger than the number of fit parameters.
 *  Assumes that the data points are independent (i.e., each point having a different x value).
 *	For power, exp by linear regression and 'Rodbard NIH Image', this is calculated for the
 *	fit actually done, not for the residuals of the original data.
 */
double CurveFitter::getFitGoodness()
{
    if (std::isnan(sumY))
        calculateSumYandY2();
    double sumMeanDiffSqr = sumY2 - sumY * sumY / sumWeights;
    double fitGoodness = 0.0;
    int degreesOfFreedom = numPoints - getNumParams();
    if (sumMeanDiffSqr > 0.0 && degreesOfFreedom > 0)
        fitGoodness = 1.0
                      - (getSumResidualsSqr() / sumMeanDiffSqr) * numPoints
                            / (double) degreesOfFreedom;

    return fitGoodness;
}

/**
 * @brief Set the maximum error. by which the sum of residuals may deviate from the true value
 *	(relative w.r.t. full sum of rediduals). Possible range: 0.1 ... 10^-16 */
void CurveFitter::setMaxError(double maxRelError)
{
    if (std::isnan(maxRelError))
        return;
    if (maxRelError > 0.1)
        maxRelError = 0.1;
    if (maxRelError < 1e-16)
        maxRelError = 1e-16; // can't be less than numerical accuracy
    this->maxRelError = maxRelError;
}

double CurveFitter::userFunction(Array &params, double dummy)
{
    double sumResidualsSqr = 0.0;
    if (numRegressionParams == 0)
    { // simply calculate sum of residuals
        for (int i = 0; i < numPoints; i++)
        {
            double fValue = f(params, xData.constData()[i]);
            double resSqr = sqr(fValue - yData.constData()[i]);
            if (!weights.isEmpty())
                resSqr *= weights.constData()[i];
            sumResidualsSqr += resSqr;
        }
    }
    else
    { // handle simple linear dependencies by linear regression:

        minimizerParamsToFullParams(params, true);
        doRegression(params);
        sumResidualsSqr = fullParamsToMinimizerParams(params);
    }
    return sumResidualsSqr;
}

/**
 * @brief For fits where linear regression is used to reduce the number of parameters handled
 *	by the Minimizer, convert Minimizer parameters to the complete set of parameters.
 *	When not for calculating regression, we use the sum of squared residuals,
 *	offset and multiplication factor stored in the extra array elements:
 *	The Minimizer stores the sum of squared residuals directly behind its last parameter.
 *	The next element is the value of the offsetParam (if any).
 *	The final element is the value of the factorParam (if any).
 */
void CurveFitter::minimizerParamsToFullParams(Array &params, bool forRegression)
{
    bool shouldTransformToSmallerParams = false;
    double offset = 0;
    double factor
        = hasSlopeParam
              ? 0
              : 1; // for regression, calculate function value without x*slope, but multiplied with 1
    double sumResidualsSqr = 0;
    if (!forRegression)
    {              // recover regression-calculated parameters from extra array elements
        int i = params.size() - 1;
        if (factorParam >= 0)
            factor = params[i--];
        if (offsetParam >= 0)
            offset = params[i];
        sumResidualsSqr
            = params[numParams
                     - numRegressionParams]; // sum of squared residuals has been calculated already
        params[numParams] = sumResidualsSqr; // ... and is now stored in its new (proper) place
    }
    // move array elements to position in array with full number of parameters
    for (int i = numParams - 1, iM = numParams - numRegressionParams - 1; i >= 0; i--)
    {
        if (i == offsetParam)
            params[i] = offset;
        else if (i == factorParam)
            params[i] = factor;
        else if (iM >= 0)
            params[i] = params[iM--];
        else
            params[i] = qNaN;
    }
    params[numParams] = sumResidualsSqr;
}

/**
 * @brief Determine sum of squared residuals with linear regression.
 *	The sum of squared residuals is written to the array element with index 'numParams',
 *	the offset and factor params (if any) are written to their proper positions in the
 *	params array */
void CurveFitter::doRegression(Array &params)
{
    double sumX = 0, sumX2 = 0, sumXY = 0; // sums for regression; here 'x' are function values
    double sumY = 0,
           sumY2 = 0; // only calculated for 'slope', otherwise we use the values calculated already
    double sumWeights = 0;
    for (int i = 0; i < numPoints; i++)
    {
        double fValue = fitType == STRAIGHT_LINE
                            ? 0
                            : f(params, xData.constData()[i]); // function value
        if (std::isnan(fValue))
        {           // check for NaN now; later we need NaN checking for division-by-zero check.
            params[numParams] = qNaN;
            return; // sum of squared residuals is NaN if any value is NaN
        }
        double w = weights.isEmpty() ? 1 : weights.constData()[i];
        sumWeights += w;

        if (hasSlopeParam)
        { // fit y = offset + slope*x + function(of other params)
            double x = xData.constData()[i];
            double y = yData.constData()[i] - fValue;
            sumX += x * w;
            sumX2 += x * x * w;
            sumXY += x * y * w;
            sumY2 += y * y * w;
            sumY += y * w;
        }
        else
        { // fit y = offset + factor * function(of other params)
            double x = fValue;
            double y = yData.constData()[i];
            sumX += fValue * w;
            sumX2 += fValue * fValue * w;
            sumXY += fValue * yData.constData()[i] * w;
        }
    }
    if (!hasSlopeParam)
    {
        sumY = this->sumY;
        sumY2 = this->sumY2;
    }
    double factor = 0;  // factor or slope
    double sumResidualsSqr = 0;
    if (offsetParam < 0)
    {                   // only 'factor' parameter, no offset:
        factor = sumXY / sumX2;
        if (std::isnan(factor) || std::isinf(factor))
            factor = 0; // all 'x' values are 0, any factor (slope) will fit
        sumResidualsSqr = sumY2 + factor * (factor * sumX2 - 2 * sumXY);
        if (sumResidualsSqr < 2e-15 * sumY2)
            sumResidualsSqr = 2e-15 * sumY2;
    }
    else
    { // full linear regression or offset only. Slope is named 'factor' here
        if (factorParam >= 0)
        {
            factor = (sumXY - sumX * sumY / sumWeights) / (sumX2 - sumX * sumX / sumWeights);
            if (restrictPower & factor <= 0) // power-law fit with (0,0) point: power must be >0
                factor = 1e-100;
            else if (std::isnan(factor) || std::isinf(factor))
                factor = 0; // all 'x' values are equal, any factor (slope) will fit
        }
        double offset = (sumY - factor * sumX) / sumWeights;
        params[offsetParam] = offset;
        sumResidualsSqr = sqr(factor) * sumX2 + sumWeights * sqr(offset) + sumY2
                          + 2 * factor * offset * sumX - 2 * factor * sumXY - 2 * offset * sumY;
        // check for accuracy problem: large difference of small numbers?
        // Don't report unrealistic or even negative values, otherwise minimization could lead
        // into parameters where we have a numeric problem
        if (sumResidualsSqr < 2e-15 * (sqr(factor) * sumX2 + sumWeights * sqr(offset) + sumY2))
            sumResidualsSqr = 2e-15 * (sqr(factor) * sumX2 + sumWeights * sqr(offset) + sumY2);
    }
    params[numParams] = sumResidualsSqr;
    if (factorParam >= 0)
        params[factorParam] = factor;
}

/**
 * @brief Convert full set of parameters to minimizer parameters and returns the sum of residuals
 * squared. The last array elements, not used by the minimizer, are the offset and factor parameters
 * (if any)
 */
double CurveFitter::fullParamsToMinimizerParams(Array &params)
{
    double offset = offsetParam >= 0 ? params[offsetParam] : 0;
    double factor = factorParam >= 0 ? params[factorParam] : 0;
    double sumResidualsSqr = params[numParams];

    for (int i = 0, iNew = 0; i < numParams;
         i++) // leave only the parameters for the minimizer in the beginning of the array
        if (i != factorParam && i != offsetParam)
            params[iNew++] = params[i];
    int i = params.size() - 1;
    if (factorParam >= 0)
        params[i--] = factor;
    if (offsetParam >= 0)
        params[i--] = offset;
    params[i--] = sumResidualsSqr;
    return sumResidualsSqr;
}

/**
 * @brief In case one or two parameters are calculated by regression and not by the minimizer:
 *	Make modified initialParams and initialParamVariations for the Minimizer
 */
void CurveFitter::modifyInitialParamsAndVariations()
{
    minimizerInitialParams = initialParams;
    minimizerInitialParamVariations = initialParamVariations;
    if (numRegressionParams
        > 0) // convert to shorter arrays with only the parameters used by the minimizer
        for (int i = 0, iNew = 0; i < numParams; i++)
            if (i != factorParam && i != offsetParam)
            {
                minimizerInitialParams[iNew] = minimizerInitialParams[i];
                minimizerInitialParamVariations[iNew] = minimizerInitialParamVariations[i];
                iNew++;
            }
}

/**
 * @brief Estimate values for initial parameters and their typical range for the built-in
 *	function types.	 For fits with modifications for regression, 'fitType' is still
 *	the type of the original (unmodified) fit.
 *	Also checks for x values being non-negative for fit types that require this,
 *	and returns false if the data cannot be fitted because of negative x.
 *	In such a case, 'errorString' contains a message about the problem. */
bool CurveFitter::makeInitialParamsAndVariations(int fitType)
{
    if (numPoints == 0)
    {
        errorString = "No data to fit";
        return false;
    }
    bool hasInitialParams = !initialParams.isEmpty();
    bool hasInitialParamVariations = !initialParamVariations.isEmpty();
    if (!hasInitialParams)
    {
        initialParams = Array(numParams);
        if (fitType == CUSTOM)
        {
            for (int i = 0; i < numParams; i++)
                initialParams[i] = 1.0;
        }
    }
    if (!hasInitialParamVariations)
        initialParamVariations = Array(numParams);
    if (fitType == CUSTOM)
    {
        for (int i = 0; i < numParams; i++)
        {
            initialParamVariations[i] = 0.1 * initialParams[i];
            if (initialParamVariations[i] == 0)
                initialParamVariations[i] = 0.01; // should not be zero
        }
        return true; // can't guess the initial parameters or initialParamVariations from the data
    }

    // Calculate some things that might be useful for predicting parameters
    double firstx = xData[0];
    double firsty = yData[0];
    double lastx = xData[numPoints - 1];
    double lasty = yData[numPoints - 1];
    double xMin = firstx, xMax = firstx; // with this initialization, the loop starts at 1, not 0
    double yMin = firsty, yMax = firsty;
    double xMean = firstx, yMean = firsty;
    double xOfMax = firstx;
    for (int i = 1; i < numPoints; i++)
    {
        xMean += xData.constData()[i];
        yMean += yData.constData()[i];
        if (xData.constData()[i] > xMax)
            xMax = xData.constData()[i];
        if (xData.constData()[i] < xMin)
            xMin = xData.constData()[i];
        if (yData.constData()[i] > yMax)
        {
            yMax = yData.constData()[i];
            xOfMax = xData.constData()[i];
        }
        if (yData.constData()[i] < yMin)
            yMin = yData.constData()[i];
    }
    xMean /= numPoints;
    yMean /= numPoints;

    double slope = (lasty - firsty) / (lastx - firstx);
    if (std::isnan(slope) || std::isinf(slope))
        slope = 0;
    double yIntercept = yMean - slope * xMean;

    // We cannot fit the following cases because we would always get NaN function values
    if (xMin < 0 && (fitType == POWER || fitType == CHAPMAN))
    {
        errorString = "Cannot fit " + fitList()[fitType] + " when x<0";
        return false;
    }
    else if (xMin < 0 && xMax > 0 && fitType == RODBARD)
    {
        errorString = "Cannot fit " + fitList()[fitType] + " to mixture of x<0 and x>0";
        return false;
    }
    else if (xMin <= 0 && fitType == LOG)
    {
        errorString = "Cannot fit " + fitList()[fitType] + " when x<=0";
        return false;
    }

    if (!hasInitialParams)
    {
        switch (fitType)
        {
        // case POLY2: case POLY3: case POLY4: case POLY5: case POLY6: case POLY7: case POLY8: {
        //  offset&slope calculated via regression; leave the others at 0

        // also for the other cases, some initial parameters are unused; only to show them with
        // 'showSettings'
        case EXPONENTIAL: // a*exp(bx)   assuming change by factor of e between xMin & xMax
            initialParams[1] = 1.0 / (xMax - xMin + 1e-100) * signum(yMean) * signum(slope);
            initialParams[0] = yMean
                               / std::exp(initialParams[1]
                                          * xMean); // don't care, done by regression
            break;
        case EXP_WITH_OFFSET: // a*exp(-bx) + c	assuming b>0, change by factor of e between xMin & xMax
        case EXP_RECOVERY:    // a*(1-exp(-bx)) + c
            initialParams[1] = 1. / (xMax - xMin + 1e-100);
            initialParams[0] = (yMax - yMin + 1e-100) / std::exp(initialParams[1] * xMean)
                                           * signum(slope) * fitType
                                       == EXP_RECOVERY
                                   ? 1
                                   : -1;    // don't care, we will do this via regression
            initialParams[2] = 0.5 * yMean; // don't care, we will do this via regression
            break;
        case EXP_RECOVERY_NOOFFSET:         // a*(1-exp(-bx))
            initialParams[1] = 1.0 / (xMax - xMin + 1e-100) * signum(yMean) * signum(slope);
            initialParams[0] = yMean
                               / std::exp(initialParams[1]
                                          * xMean);    // don't care, done by regression
            break;
        case POWER:                                    // ax^b, assume linear for the beginning
            initialParams[0] = yMean
                               / (std::abs(xMean
                                           + 1e-100)); // don't care, we will do this via regression
            initialParams[1] = 1.0;
            break;
        case LOG:                                      // a*ln(bx), assume b=e/xMax
            initialParams[0] = yMean;                  // don't care, we will do this via regression
            initialParams[1] = std::exp(1.0) / (xMax + 1e-100);
            break;
        case LOG2:                                     // y = a+b*ln(x-c)
            initialParams[0] = yMean;                  // don't care, we will do this via regression
            initialParams[1] = (yMax - yMin + 1e-100)
                               / (xMax - xMin
                                  + 1e-100);           // don't care, we will do this via regression
            initialParams[2] = jmin(0., -xMin - 0.1 * (xMax - xMin) - 1e-100);
            break;
        case RODBARD:                                  // d+(a-d)/(1+(x/c)^b)
            initialParams[0] = firsty;                 // don't care, we will do this via regression
            initialParams[1] = 1.0;
            initialParams[2] = xMin < 0 ? xMin : xMax; // better than xMean;
            initialParams[3] = lasty;                  // don't care, we will do this via regression
            break;
        case INV_RODBARD:
        case RODBARD2:                                 // c*((x-a)/(d-x))^(1/b)
            initialParams[0] = xMin - 0.1 * (xMax - xMin);
            initialParams[1] = slope >= 0 ? 1.0 : -1.0;
            initialParams[2] = yMax;                   // don't care, we will do this via regression
            initialParams[3] = xMax + (xMax - xMin);
            break;
        case GAMMA_VARIATE:
        { // // b*(x-a)^c*exp(-(x-a)/d)
            //	First guesses based on following observations (naming the 'x' coordinate 't'):
            //	t0 [a] = time of first rise in gamma curve - so use the user specified first x value
            //	tmax = t0 + c*d where tmX is the time of the peak of the curve
            //	therefore an estimate for c and d is sqrt(tm-t0)
            //	K [a] can now be calculated from these estimates
            initialParams[0] = xMin;
            double cd = xOfMax - xMin;
            if (cd < 0.1 * (xMax - xMin))
                cd = 0.1 * (xMax - xMin);
            initialParams[2] = std::sqrt(cd);
            initialParams[3] = std::sqrt(cd);
            initialParams[1]
                = yMax
                  / (std::pow(cd, initialParams[2])
                     * std::exp(-cd
                                / initialParams[3])); // don't care, we will do this via regression
            break;
        }
        case GAUSSIAN:                                // a + (b-a)*exp(-(x-c)^2/(2d^2))
            initialParams[0] = yMin; // actually don't care, we will do this via regression
            initialParams[1] = yMax; // actually don't care, we will do this via regression
            initialParams[2] = xOfMax;
            initialParams[3] = 0.39894 * (xMax - xMin) * (yMean - yMin) / (yMax - yMin + 1e-100);
            break;
        case GAUSSIAN_NOOFFSET:        // a*exp(-(x-b)^2/(2c^2))
            initialParams[0] = yMax;   // actually don't care, we will do this via regression
            initialParams[1] = xOfMax; // actually don't care, we will do this via regression
            initialParams[2] = 0.39894 * (xMax - xMin) * yMean / (yMax + 1e-100);
            break;
        case CHAPMAN:                  // a*(1-exp(-b*x))^c
            initialParams[0] = yMax;
            initialParams[2] = 1.5;    // just assuming any reasonable value
            for (int i = 1; i < numPoints; i++) // find when we reach 50% of maximum
                if (yData.constData()[i] > 0.5 * yMax)
                {                               // approximately (1-exp(-1))^1.5 = 0.5
                    initialParams[1] = 1. / xData.constData()[i];
                    break;
                }
            if (std::isnan(initialParams[1])
                || initialParams[1]
                       > 1000. / xMax)     // in case an outlier at the beginning has fooled us
                initialParams[1] = 10. / xMax;
            break;
        case ERF:                          // a+b*erf((x-c)/d)
            initialParams[0] = 0.5
                               * (yMax
                                  + yMin); // actually don't care, we will do this via regression
            initialParams[1] = 0.5 * (yMax - yMin + 1e-100)
                               * (lasty > firsty
                                      ? 1
                                      : -1); // actually don't care, we will do this via regression
            initialParams[2] = xMin
                               + (xMax - xMin) * (lasty > firsty ? yMax - yMean : yMean - yMin)
                                     / (yMax - yMin + 1e-100);
            initialParams[3] = 0.1 * (xMax - xMin + 1e-100);
            break;
            // no case CUSTOM: here, was done above
        }
    }
    if (!hasInitialParamVariations)
    {                                     // estimate initial range for parameters
        for (int i = 0; i < numParams; i++)
            initialParamVariations[i]
                = 0.1 * initialParams[i]; // default, should be overridden if it can be zero
        switch (fitType)
        {
        case POLY2:
        case POLY3:
        case POLY4:
        case POLY5:
        case POLY6:
        case POLY7:
        case POLY8:
        {
            double xFactor = 0.5 * jmax(std::abs(xMax + xMin), (xMax - xMin));
            initialParamVariations[numParams - 1] = (yMax - yMin)
                                                    / (std::pow(0.5 * (xMax - xMin), numParams - 1)
                                                       + 1e-100);
            for (int i = numParams - 2; i >= 0; i--)
                initialParamVariations[i] = initialParamVariations[i + 1] * xFactor;
            break;
        }
        case EXPONENTIAL:     // a*exp(bx)		  a (and c) is calculated by regression
        case EXP_WITH_OFFSET: // a*exp(-bx) + c
        case EXP_RECOVERY:    // a*(1-exp(-bx)) + c
            initialParamVariations[1] = 0.1 / (xMax - xMin + 1e-100);
            break;
        // case CHAPMAN:            // a*(1-exp(-b*x))^c use default (10% of value) for b, c
        // case POWER:				// ax^b; use default for b
        // case LOG:				// a*ln(bx); use default for b
        // case LOG2:				// y = a+b*ln(x-c); use default for c
        case RODBARD: // d+(a-d)/(1+(x/c)^b); a and d calculated by regression
            initialParamVariations[2] = 0.5 * jmax((xMax - xMin), std::abs(xMean));
            initialParamVariations[3] = 0.5 * jmax(yMax - yMin, std::abs(yMax));
            break;
        case INV_RODBARD: // c*((x-a)/(d-x))^(1/b); c calculated by regression
            initialParamVariations[0] = 0.01 * jmax(xMax - xMin, std::abs(xMax));
            initialParamVariations[2] = 0.1 * jmax(yMax - yMin, std::abs(yMax));
            initialParamVariations[3] = 0.1 * jmax((xMax - xMin), std::abs(xMean));
            break;
        case GAMMA_VARIATE:
        { // // b*(x-a)^c*exp(-(x-a)/d); b calculated by regression
            //	First guesses based on following observations:
            //	t0 [b] = time of first rise in gamma curve - so use the user specified first limit
            //	tm = t0 + a*B [c*d] where tm is the time of the peak of the curve
            //	therefore an estimate for a and B is sqrt(tm-t0)
            //	K [a] can now be calculated from these estimates
            initialParamVariations[0] = 0.1 * jmax(yMax - yMin, std::abs(yMax));
            double ab = xOfMax - firstx + 0.1 * (xMax - xMin + 1e-100);
            initialParamVariations[2] = 0.1 * std::sqrt(ab);
            initialParamVariations[3] = 0.1 * std::sqrt(ab);
            break;
        }
        case GAUSSIAN: // a + (b-a)*exp(-(x-c)^2/(2d^2)); a,b calculated by regression
            initialParamVariations[2] = 0.2 * initialParams[3]; //(and default for d)
            break;
        case GAUSSIAN_NOOFFSET:                                 // a*exp(-(x-b)^2/(2c^2))
            initialParamVariations[1] = 0.2 * initialParams[2]; //(and default for c)
            break;
        case ERF:                                               // a+b*erf((x-c)/d)
            initialParamVariations[2] = 0.1 * (xMax - xMin + 1e-100);
            initialParamVariations[3] = 0.5 * initialParams[3];
            break;
        }
    }
    return true;
}

/**
 * @brief Set multiplyParams and offsetParam for built-in functions. This allows us to use linear
 *	regression, reducing the number of parameters used by the Minimizer by up to 2, and
 *	improving the speed and success rate of the minimization process */
void CurveFitter::getOffsetAndFactorParams()
{
    offsetParam = -1;
    factorParam = -1;
    hasSlopeParam = false;
    switch (fitType)
    {
    case STRAIGHT_LINE:
    case POLY2:
    case POLY3:
    case POLY4:
    case POLY5:
    case POLY6:
    case POLY7:
    case POLY8:
        offsetParam = 0;
        factorParam = 1;
        hasSlopeParam = true;
        break;
    case EXPONENTIAL:     // a*exp(bx)
        factorParam = 0;
        break;
    case EXP_WITH_OFFSET: // a*exp(-bx) + c
    case EXP_RECOVERY:    // a*(1-exp(-bx)) + c
        offsetParam = 2;
        factorParam = 0;
        break;
    case EXP_RECOVERY_NOOFFSET: // a*(1-exp(-bx))
        factorParam = 0;
        break;
    case CHAPMAN:               // a*(1-exp(-b*x))^c
        factorParam = 0;
        break;
    case POWER:                 // ax^b
        factorParam = 0;
        break;
    case LOG:                   // a*ln(bx)
        factorParam = 0;
        break;
    case LOG2:                  // y = a+b*ln(x-c)
        offsetParam = 0;
        factorParam = 1;
        break;
    case RODBARD_INTERNAL: // d+a/(1+(x/c)^b)
        offsetParam = 3;
        factorParam = 0;
        break;
    case INV_RODBARD:       // c*((x-a)/(d-x))^(1/b)
        factorParam = 2;
        break;
    case GAMMA_VARIATE:     // b*(x-a)^c*exp(-(x-a)/d)
        factorParam = 1;
        break;
    case GAUSSIAN_INTERNAL: // a + b*exp(-(x-c)^2/(2d^2))
        offsetParam = 0;
        factorParam = 1;
        break;
    case GAUSSIAN_NOOFFSET: // a*exp(-(x-b)^2/(2c^2))
        factorParam = 0;
        break;
    case ERF:               // a + b*erf((x-c)/d)
        offsetParam = 0;
        factorParam = 1;
        break;
    }
    numRegressionParams = 0;
    if (offsetParam >= 0)
        numRegressionParams++;
    if (factorParam >= 0)
        numRegressionParams++;
}

/**
 * @brief calculates the sum of y and y^2 (weighted sum if we have weights) */
void CurveFitter::calculateSumYandY2()
{
    sumY = 0.0;
    sumY2 = 0.0;
    sumWeights = 0.0;
    double w = 1.0;
    for (int i = 0; i < numPoints; i++)
    {
        double y = yData.constData()[i];
        if (!weights.isEmpty())
            w = weights.constData()[i];
        sumY += y * w;
        sumY2 += y * y * w;
        sumWeights += w;
    }
}

/**
 * @brief returns whether this a fit type that acutally fits modified data with a modified function */
bool CurveFitter::isModifiedFitType(int fitType)
{
    return fitType == POWER_REGRESSION || fitType == EXP_REGRESSION || fitType == RODBARD
           || fitType == RODBARD2 || fitType == GAUSSIAN;
}

/**
 * @brief For fits don't use the original data, prepare modified data and fit type.
 *	Returns false if the data are incompatible with the fit type
 *	In that case, 'errorString' is set to a message explaining the problem
 */
bool CurveFitter::prepareModifiedFitType(int fitType)
{
    if (fitType == GAUSSIAN)
    {
        this->fitType = GAUSSIAN_INTERNAL; // different definition of parameters for using regression
        return true;
    }
    else if (fitType == RODBARD)
    {
        this->fitType = RODBARD_INTERNAL; // different definition of parameters for using regression
        return true;
    }
    else if (fitType == POWER_REGRESSION || fitType == EXP_REGRESSION)
    {
        if (fitType == POWER_REGRESSION)
        {
            xDataSave = xData;
            xData = Array(numPoints);
        }
        yDataSave = yData;
        yData = Array(numPoints);
        ySign = 0;
        numPoints = 0; // we may have lower number of points if there is a (0,0) point that we ignore
        for (int i = 0; i < xData.size(); i++)
        {
            double y = yDataSave[i];
            if (fitType == POWER_REGRESSION)
            {
                double x = xDataSave[i];
                if (x == 0 && y == 0)
                {
                    restrictPower = true;
                    continue; // ignore (0,0) point in power law
                }
                if (x <= 0)
                {
                    errorString = "Cannot fit x<=0";
                    return false;
                }
                xData[numPoints] = std::log(x);
            }
            if (ySign == 0)
                ySign = signum(y); // if unknown, determine whether y data are positive or negative
            if (y * ySign <= 0)
            {
                errorString = "Cannot fit y=0 or mixture of y>0, y<0";
                return false;
            }
            yData[numPoints] = std::log(y * ySign);
            numPoints++;
        }
        this->fitType = STRAIGHT_LINE;
    }
    else if (fitType == RODBARD2)
    { // 'Rodbard (NIH Image)' fit is inverse Rodbard function, by fitting x(y); for compatibility
      // with NIH Image
        xDataSave = xData;
        yDataSave = yData;
        xData = yDataSave; // swap
        yData = xDataSave;
        this->fitType = RODBARD_INTERNAL;
    }
    return true;
}

/**
 * @brief Get correct params and revert xData, yData if fit has been done via another function */
void CurveFitter::postProcessModifiedFitType(int fitType)
{
    if (fitType == POWER_REGRESSION
        || fitType == EXP_REGRESSION)         // ln y = ln (a*x^b) = ln a + b ln x
        finalParams[0] = ySign
                         * std::exp(
                             finalParams[0]); // or: ln (+,-)y = ln ((+,-)a*exp(bx)) = ln (+,-)a + bx
    if (fitType == GAUSSIAN)                  // a + b exp(-...) to  a + (b-a)*exp(-...)
        finalParams[1] += finalParams[0];
    else if (fitType == RODBARD || fitType == RODBARD2) // d+a/(1+(x/c)^b) to d+(a-d)/(1+(x/c)^b)
        finalParams[0] += finalParams[3];

    if (!xDataSave.isEmpty())
    {
        xData = xDataSave;
        numPoints = xData.size();
    }
    if (!yDataSave.isEmpty())
        yData = yDataSave;
    this->fitType = fitType;
}

/**
 * @brief Gets index of highest value in an array.
 *
 * @param 		   array the array.
 * @return			   Index of highest value.
 */
int CurveFitter::getMax(const Array &array)
{
    double max = array[0];
    int index = 0;
    for (int i = 1; i < array.size(); i++)
    {
        if (max < array[i])
        {
            max = array[i];
            index = i;
        }
    }
    return index;
}

/**
 * @brief 生成 ImageJ 报告中的数值统计和参数列表。
 * @return 包含公式、状态、迭代次数、SSE、标准差和决定系数的文本。
 */
QString CurveFitter::getResultString()
{
    QString report = tr("\nFormula: %1\nMacro code: %2\nStatus: %3")
                         .arg(getFormula(), getMacroCode(), getStatusString());
    if (getStatus() == Minimizer::INITIALIZATION_FAILURE)
        return report;
    if (!linearRegressionUsed)
        report += tr("\nNumber of completed minimizations: %1")
                      .arg(minimizer.getCompletedMinimizations());
    report += tr("\nNumber of iterations: %1").arg(getIterations());
    if (!linearRegressionUsed)
        report += tr(" (max: %1)").arg(minimizer.getMaxIterations());
    report += tr("\nTime: %1 ms\nSum of residuals squared: %2\nStandard deviation: %3\nR^2: "
                 "%4\nParameters:")
                  .arg(time)
                  .arg(getSumResidualsSqr(), 0, 'g', 10)
                  .arg(getSD(), 0, 'g', 10)
                  .arg(getRSquared(), 0, 'g', 10);
    const auto p = getParams();
    for (int i = 0; i < numParams; ++i)
        report += QString("\n\t%1 = %2").arg(QChar('a' + i)).arg(p[i], 0, 'g', 10);

    return report;
}

/**
 * @brief 沿用原版 Plot 的采样点数、256 灰阶特例、坐标范围和图例生成规则。
 * @param points 请求的采样点数，256 表示固定 x=0...255。
 * @return 曲线与原始散点的绘图数据，不依赖 ImageJ AWT 窗口。
 */
CurveFitter::PlotData CurveFitter::getPlot(int points)
{
    PlotData plot;
    plot.title = getFormula();
    plot.xLabel = "X";
    plot.yLabel = "Y";
    for (int i = 0; i < xData.size(); ++i)
        plot.data.append(QPointF(xData.constData()[i], yData.constData()[i]));
    if (getStatus() == Minimizer::INITIALIZATION_FAILURE || xData.isEmpty())
    {
        plot.legend = getName() + '\n' + getStatusString();
        return plot;
    }

    const int count = points == 256 ? 256 : jmin(1000, jmax(points, xData.size()));
    plot.xMin = points == 256 ? 0 : *std::min_element(xData.cbegin(), xData.cend());
    plot.xMax = points == 256 ? 255 : *std::max_element(xData.cbegin(), xData.cend());
    double low = *std::min_element(yData.cbegin(), yData.cend());
    double high = *std::max_element(yData.cbegin(), yData.cend());
    double curveLow = std::numeric_limits<double>::infinity(), curveHigh = -curveLow;
    double sampleX = plot.xMin;
    const double step = (plot.xMax - plot.xMin) / (count - 1);
    const Array params = getParams();
    for (int i = 0; i < count; ++i)
    {
        // 原版先将采样横坐标和纵坐标存为 float，再求显示范围。
        const float x = static_cast<float>(sampleX);
        sampleX += step;
        const float y = static_cast<float>(f(params, x));
        plot.curve.append(QPointF(x, y));
        if (!std::isnan(y))
        {
            curveLow = jmin(curveLow, y);
            curveHigh = jmax(curveHigh, y);
        }
    }
    const double range = high - low;
    plot.yMin = jmax(low - range, jmin(low, curveLow));
    plot.yMax = jmin(high + range, jmax(high, curveHigh));
    plot.legend = getName() + '\n' + getFormula();
    for (int i = 0; i < numParams; ++i)
        plot.legend += QString("\n%1 = %2").arg(QChar('a' + i)).arg(params[i], 0, 'g', 10);
    plot.legend += tr("\nR^2 = %1").arg(getRSquared(), 0, 'g', 5);

    return plot;
}
