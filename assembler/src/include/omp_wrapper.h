#ifndef __OMP_WRAPPER_H__
#define __OMP_WRAPPER_H__

#ifdef _OPENMP
# include <omp.h>
#else
/* Sanity check. Allow lack of OpenMP on Darwin only for now */
# ifndef __APPLE__
#  error "SPAdes requires OpenMP on this platform"
# endif
/* Provide single-threaded stubs */
# define omp_set_num_threads(x)  ((void)x)
# define omp_get_max_threads()   1
# define omp_get_thread_num()    0
# define omp_get_num_threads()   1
#endif

#endif /* __OMP_WRAPPER_H__ */
