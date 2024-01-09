#ifndef ENUMS_H
#define ENUMS_H

// TODO: capitalise for consistency
enum class designs { NS, SS, SSR };
enum class methods { mean, delta };
enum class dists { rgamma, rbeta, rnbinom, rpois, identity, rlnorm };

enum class Results {
  zero_pre, // Failure due to zero mean pre-treatment
  few_screen, // Failure due to insufficient screening positive
  few_pre, // Failure due to insufficient pre-treatment positive
  efficacy_below, // Where no test is applied but the efficacy is below the lower target
  efficacy_above, // Where no test is applied but the efficacy is greater than or equal to the lower target
  class_fail, // Unable to apply the test/classification due to 100% reduction or something similar
  resistant, // Where a test is applied
  low_resistant, // Where a test is applied
  inconclusive, // Where a test is applied
  susceptible // Where a test is applied
};

// TODO: unit test over all Results
constexpr const char* ResultToString(Results result)
{
    switch (result)
    {
        case Results::zero_pre: return "FailZeroPre";
        case Results::few_screen: return "FailPositiveScreen";
        case Results::few_pre: return "FailPositivePre";
        case Results::efficacy_below: return "EfficacyBelowLT"; // Only for mean method
        case Results::efficacy_above: return "EfficacyAboveLT"; // Only for mean method (actually above or equal to)
        case Results::class_fail: return "ClassifyFail"; // This and below only for delta method
        case Results::resistant: return "Resistant";
        case Results::low_resistant: return "LowResistant";
        case Results::inconclusive: return "Inconclusive";
        case Results::susceptible: return "Susceptible";
        default: Rcpp::stop("Unimplemented result in ResultToString");
    }
}

constexpr const bool ResultIsSuccess(Results result)
{
    switch (result)
    {
        case Results::zero_pre: return false;
        case Results::few_screen: return false;
        case Results::few_pre: return false;
        case Results::efficacy_below: return true; // Only for mean method
        case Results::efficacy_above: return true; // Only for mean method (actually above or equal to)
        case Results::class_fail: return true; // NOTE: this is a success in that means can be calculated
        case Results::resistant: return true;
        case Results::low_resistant: return true;
        case Results::inconclusive: return true;
        case Results::susceptible: return true;
        default: Rcpp::stop("Unimplemented result in ResultIsSuccess");
    }
}

// TODO: unit test against above
constexpr const bool ResultIsSuccess(const int result)
{
    switch (result)
    {
        case 0L: return false; // Results::zero_pre
        case 1L: return false; // Results::few_screen
        case 2L: return false; // Results::few_pre
        case 3L: return true; // Results::efficacy_below
        case 4L: return true; // Results::efficacy_above
        case 5L: return true; // Results::class_fail
        case 6L: return true; // Results::resistant
        case 7L: return true; // Results::low_resistant
        case 8L: return true; // Results::inconclusive
        case 9L: return true; // Results::susceptible
        default: Rcpp::stop("Unimplemented result in ResultIsSuccess");
    }
}


// Note: constexpr std::array<const char*, 9L> works on arm64 but not x86_64 ?!?
// Note: std::string can't be const on Windows :(
/*
static const std::vector<std::string> GetResultsLevels()
{
  const std::vector<std::string> ResultsLevels = {
      ResultToString(Results::zero_pre),
      ResultToString(Results::few_screen),
      ResultToString(Results::few_pre),
      ResultToString(Results::efficacy_below),
      ResultToString(Results::efficacy_above),
      ResultToString(Results::class_fail),
      ResultToString(Results::resistant),
      ResultToString(Results::low_resistant),
      ResultToString(Results::inconclusive),
      ResultToString(Results::susceptible)
    };
    return ResultsLevels;
}
*/


// For testing costs with fixed data:
#define TESTING() constexpr int s_testing = 0L; constexpr bool t_testing = false

// For passing around:
struct CountParams
{
  double intercept;
  double coefficient;
  double add;
  double mult;
  int min_pos_screen;
  int min_pos_pre;
  double tail;
  double Teff;
  double Tlow;
  double inclusion_prob = 1.0;
  double retention_prob_screen = 1.0;
  double retention_prob_pre = 1.0;
};

template <typename T>
T coin_toss(const double probability)
{
  const double res = R::rbinom(1.0, probability);
  return static_cast<T>(res);
}


#endif // ENUMS_H
