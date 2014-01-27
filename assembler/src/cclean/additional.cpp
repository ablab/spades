#ifndef ADDITIONAL_CPP
#define ADDITIONAL_CPP
namespace additional {

  enum WorkModeType {
    NONE = 0,
    SIGNLE_END = 1,
    SINGLE_END_Q = 2,
    BRUTE_SIMPLE = 3,
    BRUTE_WITH_Q = 4
  };

  constexpr double MatchScore = 0.6;
  constexpr double MismatchScore = 100;

  // end of namespace additional
}
#endif // ADDITIONAL_CPP
