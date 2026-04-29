/** \brief Record running statistics on a logarithmic scale.
 *
 * This class is primarily meant to record statistics on the time steps
 * `dt` and magnetization variations per time step `maxdu`. As these
 * quantities may vary over multiple orders of magnitude, the statistics
 * of the logarithms are likely to be more relevant that the statistics
 * on the quantities themselves.
 *
 * For the sake of numerical stability, this class uses [Welford's
 * online algorithm][1].
 *
 * [1]:
 * https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford%27s_online_algorithm
 */

#include <cmath>

class LogStats
    {
public:
    /** Add a sample to the statistics. */
    void add(double x)
        {
        x = std::log(x);
        n += 1;
        double delta1 = x - m;
        m += delta1 / n;
        double delta2 = x - m;
        s += delta1 * delta2;
        }

    /** Return the count of samples added so far. */
    long count() const { return n; }

    /**
     * Return the geometric mean, i.e. the exponential of the mean of
     * the logarithms.
     *
     * For a distribution that is wide on a logarithmic scale, this is
     * usually a better indication of what a “typical” value is. For a
     * log-normal distribution, this is the same as the median.
     */
    double mean() const { return std::exp(m); }

    /**
     * Return the sample standard deviation of the logarithm.
     *
     * For a narrow distribution, this can be interpreted as a “relative
     * width”. For example, the value 0.05 means a width roughly equal
     * to 5% of the mean. Larger values can be interpreted as expressing
     * the width via multiplicative factors. For example, the value 2.3
     * (i.e. `log(10)`) can be interpreted as a distribution where most
     * samples are within a factor 10 of the geometric mean.
     */
    double stddev() const { return std::sqrt(s / n); }

private:
    long n = 0;   /**< sample count */
    double m = 0; /**< sample mean */
    double s = 0; /**< sum of squares of deviations */
    };
