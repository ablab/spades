#ifndef ADDITIONAL_CPP
#define ADDITIONAL_CPP
namespace additional {

  enum WorkModeType {
    NONE = 0,
    SIMPLE = 1,
    BRUTE_SIMPLE = 3,
    BRUTE_DEEP = 4
  };

  const double BruteMatchScore = 0.6;
  const double BruteMismatchScore = 10;

  // end of namespace additional
}
#endif // ADDITIONAL_CPP
