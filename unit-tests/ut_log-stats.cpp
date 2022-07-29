/*
 * Test the class LogStats defined in log-stats.h.
 */

#define BOOST_TEST_MODULE logStatsTest

#include <boost/test/unit_test.hpp>
#include "ut_config.h"
#include "log-stats.h"
#include <list>
#include <numeric>
#include <random>

/* Naive implementation to compare against. */
class NaiveLogStats : private std::list<double>
{
public:
    void add(double x) { push_back(std::log(x)); }
    long count() const { return size(); }
    double mean() const { return std::exp(meanLog()); }
    double stddev() const {
        double m = meanLog();
        double sum = std::accumulate(begin(), end(), 0.0,
            [m](double s, double x){ return s + (x - m) * (x - m); }
        );
        return std::sqrt(sum / size());
    }
private:
    double meanLog() const {
        return std::accumulate(begin(), end(), 0.0) / size();
    }
};

BOOST_AUTO_TEST_SUITE(ut_log_stats)

BOOST_AUTO_TEST_CASE(log_stats, * boost::unit_test::tolerance(1e-12))
{
    LogStats stats;
    NaiveLogStats reference;
    std::mt19937 generator(my_seed());
    std::normal_distribution<> normal(std::log(0.1), 3);

    for (int i = 0; i < 10000; ++i) {
        double x = std::exp(normal(generator));
        stats.add(x);
        reference.add(x);
    }

    BOOST_TEST(stats.count() == reference.count());
    BOOST_TEST(stats.mean() == reference.mean());
    BOOST_TEST(stats.stddev() == reference.stddev());
}

BOOST_AUTO_TEST_SUITE_END()
