#include "minimizer.h"
#include "numericcompat.h"

#include <exception>
#include <stdexcept>
#include <thread>

using namespace ImageJNumeric;

////////////////////////////////////////////////////////////////////////
// class Minimizer

Minimizer::JavaRandom::JavaRandom(int value)
    : seed((static_cast<std::uint64_t>(static_cast<std::int64_t>(value)) ^ 0x5DEECE66DULL)
           & ((1ULL << 48) - 1))
{
}

std::uint32_t Minimizer::JavaRandom::next(int bits)
{
    seed = (seed * 0x5DEECE66DULL + 0xBULL) & ((1ULL << 48) - 1);
    return static_cast<std::uint32_t>(seed >> (48 - bits));
}

float Minimizer::JavaRandom::nextFloat()
{
    return next(24) / static_cast<float>(1 << 24);
}

double Minimizer::JavaRandom::nextDouble()
{
    const auto high = static_cast<std::uint64_t>(next(26)) << 27;
    return (high + next(27)) / static_cast<double>(1ULL << 53);
}

QStringList Minimizer::statusStrings()
{
    return {tr("Success"),
            tr("Initialization failure; no result"),
            tr("Aborted"),
            tr("Re-initialization failure (inaccurate result?)"),
            tr("Max. no. of iterations reached (inaccurate result?)"),
            tr("Max. no. of restarts reached (inaccurate result?)")};
}

void Minimizer::setFunction(UserFunction function, int count)
{
    if (!function || count < 0)
        throw std::invalid_argument("Invalid minimization function or dimension");
    if (maxIter <= 0)
    {
        maxIter = ITER_FACTOR * count * count;
        if (maxRestarts > 0)
            maxIter *= 2;
    }
    userFunction = std::move(function);
    numParams = count;
    numVertices = count + 1;
}

void Minimizer::setMaxError(double relative, double absolute)
{
    maxRelError = relative;
    maxAbsError = absolute;
}

void Minimizer::setStatusAndEsc(const QString &prefix, bool escape)
{
    ijStatusString = prefix;
    checkEscape = escape;
}

void Minimizer::setStatusCallbacks(StatusCallback display, EscapeCallback escape)
{
    statusCallback = std::move(display);
    escapeCallback = std::move(escape);
}

void Minimizer::reportProgress()
{
    if (ijStatusString.isNull() && !checkEscape)
        return;
    std::lock_guard<std::mutex> lock(progressMutex);
    const int iterations = totalNumIter;
    if (iterations <= nextIterationForStatus)
        return;
    const auto elapsed = timer.elapsed();
    nextIterationForStatus = iterations + static_cast<int>(iterations * 500LL / (elapsed + 1));
    if (elapsed <= 1000)
        return;
    if (checkEscape && escapeCallback && escapeCallback())
    {
        status = ABORTED;
        if (statusCallback)
            statusCallback(ijStatusString + tr(" ABORTED"));
        return;
    }
    if (!ijStatusString.isNull() && statusCallback)
        statusCallback(ijStatusString + QString::number(iterations) + tr(" (%1 max)").arg(maxIter)
                       + (checkEscape ? tr(" ESC to stop") : QString()));
}

int Minimizer::minimize(const Array &initialParams, const Array &variations)
{
    status = SUCCESS;
    resultsVector.clear();
    collectResults = true;
    timer.start();
    int loops = maxRestarts + 1;
    if (useSingleThread)
        loops *= 2;
    for (int i = 0; i < loops; ++i)
    {
        std::exception_ptr workerException;
        std::thread second;
        if (maxRestarts > 0 && !useSingleThread)
            second = std::thread([&, i] {
                try
                {
                    minimizeOnce(initialParams, variations, randomSeed + 1000000 + i);
                }
                catch (...)
                {
                    workerException = std::current_exception();
                    status = ABORTED;
                }
            });
        try
        {
            minimizeOnce(initialParams, variations, randomSeed + i);
        }
        catch (...)
        {
            status = ABORTED;
            if (second.joinable())
                second.join();
            throw;
        }
        if (second.joinable())
            second.join();
        if (workerException)
            std::rethrow_exception(workerException);
        if (resultsVector.isEmpty() && result.isEmpty())
            return status;
        if (result.isEmpty())
            result = resultsVector[0];
        for (auto &r : resultsVector)
            if (value(r) < value(result))
                result = r;
        if (status != SUCCESS && status != REINITIALIZATION_FAILURE
            && status != MAX_ITERATIONS_EXCEEDED)
            return status;
        if (totalNumIter >= maxIter)
            return MAX_ITERATIONS_EXCEEDED;
        for (int ir = 0; ir < resultsVector.size(); ++ir)
            if (!belowErrorLimit(value(resultsVector[ir]), value(result), 1.0))
            {
                resultsVector.removeAt(ir);
                --ir;
            }
        if (resultsVector.size() >= 2)
            return SUCCESS;
    }
    if (!ijStatusString.isNull() && statusCallback)
        statusCallback(QString());
    return maxRestarts > 0 ? MAX_RESTARTS_EXCEEDED : status.load();
}

int Minimizer::minimizeOnce(const Array &initialParams, const Array &variations)
{
    status = SUCCESS;
    timer.start();
    minimizeOnce(initialParams, variations, randomSeed);
    return status;
}

Minimizer::Array Minimizer::getParams()
{
    if (result.isEmpty())
        result = Array(numParams + 1 + numExtraArrayElements, qNaN);
    return result;
}

double Minimizer::getFunctionValue()
{
    if (result.isEmpty())
        result = Array(numParams + 1, qNaN);
    return value(result);
}

void Minimizer::minimizeOnce(const Array &initialParams, const Array &variations, int seed)
{
    JavaRandom random(seed);
    Matrix simp = makeSimplex(initialParams, variations, random);
    if (simp.isEmpty())
    {
        status = wasInitialized ? REINITIALIZATION_FAILURE : INITIALIZATION_FAILURE;
        return;
    }
    wasInitialized = true;
    int best = minimizeSimplex(simp);
    double bestValue = value(simp[best]);
    bool failed = false;
    while (status == SUCCESS || status == REINITIALIZATION_FAILURE)
    {
        Array newVariations = makeNewParamVariations(simp, best, initialParams, variations);
        if (!reInitializeSimplex(simp, best, newVariations, random))
        {
            failed = true;
            break;
        }
        best = minimizeSimplex(simp);
        if (belowErrorLimit(value(simp[best]), bestValue, 2.0))
            break;
        bestValue = value(simp[best]);
    }
    if (failed)
        status = REINITIALIZATION_FAILURE;
    else if (status == SUCCESS || status == REINITIALIZATION_FAILURE)
        ++numCompletedMinimizations;
    std::lock_guard<std::mutex> lock(resultsMutex);
    if (collectResults)
        resultsVector.append(simp[best]);
    else
        result = simp[best];
}

int Minimizer::minimizeSimplex(Matrix &simp)
{
    QVector<int> worstNextBestArray = QVector<int>(
        3); // used to pass these indices from 'order' function
    Array center = Array(numParams + 1
                         + numExtraArrayElements); // center of all vertices except worst
    Array reflected = Array(numParams + 1 + numExtraArrayElements); // the 1st new vertex to try
    Array secondTry = Array(numParams + 1 + numExtraArrayElements); // expanded or contracted vertex

    order(simp, worstNextBestArray);
    int worst = worstNextBestArray[WORST];
    int nextWorst = worstNextBestArray[NEXT_WORST];
    int best = worstNextBestArray[BEST];

    int thisNumIter = 0;
    while (true)
    {
        totalNumIter++; // global count over all threads
        thisNumIter++;  // local count for this minimize call
        // THE MINIMIZAION ALGORITHM IS HERE
        do
        {
            getCenter(simp, worst, center); // centroid of vertices except worst
            // Reflect worst vertex through centriod of not-worst
            getVertexAndEvaluate(center, simp[worst], -C_REFLECTION, reflected);
            if (value(reflected) <= value(simp[best]))
            { // if it's better than the best...
                // try expanding it
                getVertexAndEvaluate(center, simp[worst], -C_EXPANSION, secondTry);
                if (value(secondTry) <= value(reflected))
                {
                    copyVertex(secondTry,
                               simp[worst]); // if expanded is better than reflected, keep it

                    break;
                }
            }
            if (value(reflected) < value(simp[nextWorst]))
            {
                copyVertex(reflected, simp[worst]); // keep reflected if better than 2nd worst

                break;
            }
            else if (value(reflected) < value(simp[worst]))
            {
                // try outer contraction
                getVertexAndEvaluate(center, simp[worst], -C_CONTRACTION, secondTry);
                if (value(secondTry) <= value(reflected))
                {
                    copyVertex(secondTry, simp[worst]); // keep outer contraction

                    break;
                }
            }
            else if (value(reflected) > value(simp[worst]) || std::isnan(value(reflected)))
            {
                // else inner contraction
                getVertexAndEvaluate(center, simp[worst], C_CONTRACTION, secondTry);
                if (value(secondTry) < value(simp[worst]))
                {
                    copyVertex(secondTry, simp[worst]); // keep contracted if better than 2nd worst

                    break;
                }
            }
            // if everything else has failed, contract simplex in on best
            shrinkSimplexAndEvaluate(simp, best);

            break;
        }
        while (false); // one simplex step
        bool checkParamResolution = // if new 'worst' is not close to 'best', don't check any further
            !paramResolutions.isEmpty() && belowResolutionLimit(simp[worst], simp[best]);
        order(simp, worstNextBestArray);
        worst = worstNextBestArray[WORST];
        nextWorst = worstNextBestArray[NEXT_WORST];
        best = worstNextBestArray[BEST];

        if (checkParamResolution)
            if (belowResolutionLimit(
                    simp, best)) // check whether all parameters are within the resolution limit
                break;
        if (belowErrorLimit(value(simp[best]), value(simp[worst]), 4.0))
        {                        // make sure we are at the minimum:
            getCenter(simp, -1, secondTry);        // check center of the simplex
            evaluate(secondTry);
            if (value(secondTry) < value(simp[best]))
                copyVertex(secondTry, simp[best]); // better than best: keep
        }
        if (belowErrorLimit(value(simp[best]), value(simp[worst]),
                            4.0))                  // no large spread of values
            break;                                 // looks like we are at the minimum
        if (totalNumIter > maxIter || thisNumIter > 4 * (maxIter / 10))
            status = MAX_ITERATIONS_EXCEEDED;
        if (status != SUCCESS)
            break;
        reportProgress();
    }

    return best;
}

/**
 * @brief Move along the line between center and worst, where +1 corresponds to worstVertex
 *  The new point is written to'newVertex' and target function is evaluated there
 */
void Minimizer::getVertexAndEvaluate(Array &center,
                                     Array &worstVertex,
                                     double howFar,
                                     Array &newVertex)
{
    for (int i = 0; i < numParams; i++)
        newVertex[i] = (1. - howFar) * center[i] + howFar * worstVertex[i];
    evaluate(newVertex);
}

/**
 * @brief Get center of all vertices except one to exclude
 *  (may be -1 to exclude none)
 *  Does not care about function values (i.e., last array elements) */
void Minimizer::getCenter(Matrix &simp, int excludeVertex, Array &center)
{
    fill(center, 0.);
    int nV = 0;
    for (int v = 0; v < numVertices; v++)
        if (v != excludeVertex)
        {
            for (int i = 0; i < numParams; i++)
                center[i] += simp[v][i];
            nV++;
        }
    double norm = 1.0 / nV;
    for (int i = 0; i < numParams; i++)
        center[i] *= norm;
}

void Minimizer::shrinkSimplexAndEvaluate(Matrix &simp, int best)
{
    for (int v = 0; v < numVertices; v++)
        if (v != best)
        {
            for (int i = 0; i < numParams; i++)
                simp[v][i] = C_SHRINK * simp[v][i] + (1 - C_SHRINK) * simp[best][i];
            evaluate(simp[v]);
        }
}

/**
 * @brief Whether the highest and lowest values fulfill the error limit criterion.
 *  Error limits are reduced by a factor of 'sensitivity' */
bool Minimizer::belowErrorLimit(double highest, double lowest, double sensitivity)
{
    double absError = sensitivity * std::abs(highest - lowest);
    double relError = absError / (jmax(std::abs(highest), std::abs(lowest)) + 1e-100);
    return relError < maxRelError || absError < maxAbsError;
}

/**
 * @brief Initialise the simplex and evaluate it. Returns null on failure. */
Minimizer::Matrix Minimizer::makeSimplex(const Array &initialParams,
                                         const Array &initialParamVariations,
                                         JavaRandom &random)
{
    Matrix simp = Matrix(numVertices, Array(numParams + 1 + numExtraArrayElements));

    if (!initialParams.isEmpty())
    {
        arraycopy(initialParams, 0, simp[0], 0, jmin(initialParams.size(), numParams));
    }
    evaluate(simp[0]);
    if (std::isnan(value(simp[0])))
    {
        findValidInitalParams(simp[0], initialParamVariations, random);
    }
    if (std::isnan(value(simp[0])))
    {
        return {};
    }
    if (initializeSimplex(simp, initialParamVariations, random))
        return simp;
    else
    {
        return {};
    }
}

/**
 * @brief Whether the distance between the best and any other simplex vertex
 *  is less than the resolution of the parameters */
bool Minimizer::belowResolutionLimit(Matrix &simp, int best)
{
    for (int v = 0; v < numVertices; v++)
        if (v != best && !belowResolutionLimit(simp[v], simp[best]))
            return false;
    return true;
}

bool Minimizer::belowResolutionLimit(Array &vertex1, Array &vertex2)
{
    for (int i = 0; i < numParams; i++)
        if (std::abs(vertex1[i] - vertex2[i]) >= paramResolutions[i])
            return false;
    return true;
}

/**
 * @brief Find initial parameters not yielding a result of NaN in case those given yield NaN
 *  Called with params containing the initial parameters that have been tried previously */
void Minimizer::findValidInitalParams(Array &params,
                                      const Array &initialParamVariations,
                                      JavaRandom &random)
{
    int maxAttempts = 50 * numParams
                      * numParams; // max number of attempts to find params that do not lead to NaN
    double rangeFactor = 1; // will gradually become larger to handle different orders of magnitude
    double rangeMultiplyLog = std::log(1e20)
                              / (maxAttempts
                                 - 1);    // will try up to 1e-20 to 1e20*initialParamVariations
    Array firstParams = Array(numParams); // remember starting params (which may be modified)
    Array variations = Array(numParams);  // new values of parameter variations
    for (int i = 0; i < numParams; i++)
    {
        firstParams[i] = std::isnan(params[i]) ? 0 : params[i];
        variations[i] = !initialParamVariations.isEmpty() ? initialParamVariations[i]
                                                          : 0.1 * firstParams[i];
        if (std::isnan(variations[i]) || std::abs(variations[i]) < 1e-10
            || std::abs(variations[i]) > 1e10)
            variations[i] = 0.1;
    }
    for (int attempt = 0; attempt < maxAttempts; attempt++)
    {
        for (int i = 0; i < numParams; i++)
        {
            double multiplier = attempt < maxAttempts / 10
                                    ? 1
                                    : // after a while try different orders of magnitude
                                    std::exp(rangeMultiplyLog * attempt * 2
                                             * (random.nextDouble() - 0.5));
            params[i] = multiplier
                        * (firstParams[i] + 2 * (random.nextDouble() - 0.5) * variations[i]);
        }
        evaluate(params);

        if (!std::isnan(value(params)))
            return; // found a valid parameter set
    }
}

/**
 * @brief Reinitialize an existing simplex: Create new vertices around the best one, keeping the
 * rough size of the simplex. This helps to avoid premature termination or cases where the simplex
 * has become degenerate. All points of the new simplex are evaluated already. Returns false on
 * error (if it could not create a simplex not resulting in NaN). */
bool Minimizer::reInitializeSimplex(Matrix &simp,
                                    int bestVertexNumber,
                                    const Array &paramVariations,
                                    JavaRandom &random)
{
    if (bestVertexNumber != 0)
    {
        Array swap = simp[0];
        simp[0] = simp[bestVertexNumber];
        simp[bestVertexNumber] = swap;
    }
    return initializeSimplex(simp, paramVariations, random);
}

/**
 * @brief Initialize simplex vertices except the number 0, which will be taken as starting point.
 *  All points of the new simplex are evaluated already.
 *  Returns false on error (if it could not create a simplex not resulting in NaN).
 */
bool Minimizer::initializeSimplex(Matrix &simp, const Array &paramVariations, JavaRandom &random)
{
    Array variations = Array(numParams); // improved values of parameter variations
    for (int i = 0; i < numParams; i++)
    {
        double range
            = !paramVariations.isEmpty() && i < paramVariations.size()
                  ? paramVariations[i]
                  : 0.1 * std::abs(simp[0][i]); // if nothing else given, based on initial parameter
        if (range < 1e-100 || std::isnan(range))
            range = 0.01;                       // range must be nonzero
        if (std::abs(range / simp[0][i]) < 1e-10)
            range = std::abs(simp[0][i] * 1e-10); // range must be more than very last digits
        variations[i] = range;
    }
    int maxAttempts = 100
                      * numParams; // max number of attempts to find params that do not lead to NaN
    for (int v = 1; v < numVertices; v++)
    {
        int numTries = 0;
        do
        { // try finding a vertex that does not yield NaN
            if (numTries++ > maxAttempts)
                return false;
            // Create a vector orthogonal to all others in a coordinate system where every value is
            // normalized to its respective variations value. We first use this coordinate system
            // where every parameter is in the [-1, +1] range. Don't orthogonalize after getting
            // only NaN function values in many attempts; it might be easier to find viable
            // parameters without normalization.
            for (int i = 0; i < numParams; i++)
                simp[v][i] = 2 * random.nextFloat() - 1.0;
            if (numTries
                < maxAttempts
                      / 3) // (don't orthogonalize if finding valid params was very difficult

                for (int v2 = 1; v2 < v; v2++)
                {          // to avoid a degenerate simplex, make orthogonal to all others
                    double lengthSqr = 0;
                    double innerProduct = 0;
                    for (int i = 0; i < numParams; i++)
                    {
                        double x2 = (simp[v2][i] - simp[0][i])
                                    / variations[i]; // convert to [-1, +1] range
                        lengthSqr += x2 * x2;
                        innerProduct += x2 * simp[v][i];
                    }
                    for (int i = 0; i < numParams;
                         i++) // subtract component parallel to the other vector
                        simp[v][i] -= ((simp[v2][i] - simp[0][i]) / variations[i])
                                      * (innerProduct / lengthSqr);
                }
            double sumSqr = 0; // normalize the 'excursion vector' to length = 1
            for (int i = 0; i < numParams; i++)
                sumSqr += simp[v][i] * simp[v][i];
            if (sumSqr < 1e-14)
                continue;   // bad luck with random numbers, excursion vector almost zero
            if (numTries > 2 && sumSqr > 1e-6)
                sumSqr = 1; // don't normalize if attempts with normalized were unsuccessful
            for (int i = 0; i < numParams; i++)
                simp[v][i] = simp[0][i] + simp[v][i] / std::sqrt(sumSqr) * variations[i];
            evaluate(simp[v]);
        }
        while (std::isnan(value(simp[v]))
               && (status == SUCCESS || status == REINITIALIZATION_FAILURE));
    }

    return true;
}

/**
 * @brief Calculate ranges for parameter variations for re-initializing:
 *  Note that the simplex should have adapted in size to the useful range of each
 *  parameter, but in special cases it might have become degenerate [i.e. the
 *  (hyper-)volume spanned by the vertices might be zero; in such a case it may
 *  happen that one parameter has the same value for all vertices of the simplex].
 *
 *  For each parameter, as a variation value we use the typical difference w.r.t. the best
 *  vertex of the simplex.
 *  To avoid a degenerate simplex, we have to know what the the relative orders of magnitude of
 *  the parameter variations should be.  If initialParamVariations are given, we use them as
 *  indication of the typoical relative values; otherwise we use the parameter values themselves
 *  and assume that the variation of each parameter should be roughly the same fraction of the value.
 *  If the variation of a parameter in the simplex w.r.t. this estimated variation is much
 *  smaller (<10^3) than typical for the parameters, we increase the variation for this parameter.
 */
Minimizer::Array Minimizer::makeNewParamVariations(Matrix &simp,
                                                   int bestVertexNumber,
                                                   const Array &initialParams,
                                                   const Array &initialParamVariations)
{
    Array paramVariations = Array(numParams); // will be the return value
    Array relatedTo = Array(
        numParams);           // parameter variation should be related to (in size) these values
    double logTypicalRelativeVariation = 0;
    for (int i = 0; i < numParams; i++)
    {                         // for all parameters
        double variation = 0; // roughly size of simplex in direction of that parameter
        for (int v = 0; v < numVertices; v++)
            if (v != bestVertexNumber)
            {
                double delta = std::abs(simp[v][i] - simp[bestVertexNumber][i]);
                paramVariations[i] += delta * delta;
            }
        paramVariations[i]
            = 10 * std::sqrt(paramVariations[i]); // make the new simplex larger than the old one
        relatedTo[i] = !initialParamVariations.isEmpty()
                               && initialParamVariations.size() >= numParams
                           ? initialParamVariations[i]
                           : // relate to initialParamVariations (if given),
                           jmax(std::abs(initialParams.value(i)),
                                std::abs(simp[bestVertexNumber][i]));
        double logRelativeVariation = paramVariations[i] > relatedTo[i]
                                          ? 0
                                          : std::log(paramVariations[i] / relatedTo[i]);
        logTypicalRelativeVariation += logRelativeVariation;
    }
    logTypicalRelativeVariation /= numParams;
    double typicalRelativeVariation = std::exp(logTypicalRelativeVariation);
    double WORST_RATIO
        = 1e-3; // parameter variation should not be lower than typical value by more than this
    for (int i = 0; i < numParams; i++)
    {
        if (paramVariations[i] < relatedTo[i]
            && paramVariations[i] / relatedTo[i] < typicalRelativeVariation * WORST_RATIO)
            paramVariations[i] = relatedTo[i] * typicalRelativeVariation * WORST_RATIO;
        if (!paramResolutions.isEmpty() && paramVariations[i] < 10 * paramResolutions[i])
            paramVariations[i]
                = 10
                  * paramResolutions[i]; // parameter variation must be well above numeric noise limit
    }
    return paramVariations;
}

/**
 * @brief Evaluate the target function for the parameters given by 'vertex'
 *  The function value is stored into the last (extra) array element.
 */
void Minimizer::evaluate(Array &vertex)
{
    vertex[numParams] = userFunction(vertex, 0);
}

/**
 * @brief Function value of a vertex after call to 'evaluate'
 *  (stored it as last array element) */
double Minimizer::value(Array &vertex)
{
    return vertex[numParams];
}

/**
 * @brief Find the simplex vertices with the worst, nextWorst and best values
 *  of the target function */
void Minimizer::order(Matrix &simp, QVector<int> &worstNextBestArray)
{
    int worst = 0, nextWorst = 0, best = 0;
    for (int i = 0; i < numVertices; i++)
    {
        if (value(simp[i]) < value(simp[best]))
            best = i;
        if (value(simp[i]) > value(simp[worst]))
            worst = i;
    }
    nextWorst = best;
    for (int i = 0; i < numVertices; i++)
        if (i != worst && value(simp[i]) > value(simp[nextWorst]))
            nextWorst = i;
    worstNextBestArray[WORST] = worst;
    worstNextBestArray[NEXT_WORST] = nextWorst;
    worstNextBestArray[BEST] = best;
}
