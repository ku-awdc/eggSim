#ifndef ENUMS_H
#define ENUMS_H

// TODO: capitalise for consistency
enum class designs { NS, SS, SSR };
enum class methods { mean, delta };
enum class dists { rgamma, rbeta, rnbinom, rpois, identity, rlnorm };

enum class Results { 
  few_screen, // Failure due to insufficient screening positive
  few_pre, // Failure due to insufficient pre-treatment positive
  zero_pre, // Failure due to zero mean pre-treatment
  success, // Where no test is applied and the two above don't apply
  zero_post, // Unable to apply the test due to 100% reduction
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
        case Results::few_screen: return "Failure due to insufficient screening positive";
        case Results::few_pre: return "Failure due to insufficient pre-treatment positive";
        case Results::zero_pre: return "Failure due to zero mean pre-treatment";
        case Results::success: return "Successful result (no classification)";
        case Results::zero_post: return "Unable to classify due to 100% reduction";
        case Results::resistant: return "Resistant classification";
        case Results::low_resistant: return "Low Resistant classification";
        case Results::inconclusive: return "Inconclusive classification";
        case Results::susceptible: return "Susceptible classification";
        default: Rcpp::stop("Unimplemented result in ResultToString");
    }
}

// For testing costs with fixed data:
#define TESTING() constexpr int s_testing = 0L

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
};


#endif // ENUMS_H
