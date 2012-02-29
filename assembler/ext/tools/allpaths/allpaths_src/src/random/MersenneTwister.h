/* Below is a slightly modified version of the code by 
   Nishimura and Mastumoto that now compiles cleanly under g++ and 
   also works as a library [main commented out]. 
   
   The body is in MersenneTwister.cc.  This just contains the headers.
   Please see the *.cc file for all the comments.
   
   */

// Import longlong typedef.
#include "system/Types.h"

/* initializes mt[NN] with a seed */
void init_genrand64(ulonglong seed);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array64(ulonglong init_key[], ulonglong key_length);

/* generates a random number on [0, 2^64-1]-interval */
ulonglong genrand64_int64(void);

/* generates a random number on [0, 2^63-1]-interval */
longlong genrand64_int63(void);

/* generates a random number on [0,1]-real-interval */
double genrand64_real1(void);

/* generates a random number on [0,1)-real-interval */
double genrand64_real2(void);

/* generates a random number on (0,1)-real-interval */
double genrand64_real3(void);

/* generates a random Gaussian distribution using Mersenne Twister. */
double genrand64_Box_Mueller_Gaussian(double offset = 0.0, double sigma = 1.0);

