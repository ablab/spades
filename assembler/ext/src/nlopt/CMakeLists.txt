project(nlopt CXX C)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${EXT_DIR}/include)

set(NLOPT_SRCS "api/general.c"
               "api/optimize.c"
               "api/options.c"
               "api/deprecated.c"
               "auglag/auglag.c"
               "bobyqa/bobyqa.c"
               "cdirect/cdirect.c"
               "cdirect/hybrid.c"
               "cobyla/cobyla.c"
               "crs/crs.c"
               "esch/esch.c"
               "isres/isres.c"
               "mma/mma.c"
               "mma/ccsa_quadratic.c"
               "mlsl/mlsl.c"
               "neldermead/nldrmd.c"
               "neldermead/sbplx.c"
               "newuoa/newuoa.c"
               "praxis/praxis.c"
               "util/redblack.c"
               "util/rescale.c"
               "util/stop.c"
               "util/timer.c"
               "util/qsort_r.c"
               "util/sobolseq.c"
               "util/mt19937ar.c"
)

add_library(nlopt STATIC ${NLOPT_SRCS})
