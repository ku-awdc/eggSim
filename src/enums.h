#ifndef ENUMS_H
#define ENUMS_H

enum class designs { NS, SS, SSR };
enum class methods { custom, delta };
enum class dists { rgamma, rbeta, rnbinom, rpois, identity, rlnorm };

// For testing costs with fixed data:
#define TESTING() constexpr int s_testing = 0L

#endif // ENUMS_H
