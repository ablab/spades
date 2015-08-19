//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __OMP_WRAPPER_H__
#define __OMP_WRAPPER_H__

#ifdef _OPENMP
# include <omp.h>
#else
/* Provide single-threaded stubs */
# define omp_set_num_threads(x)  ((void)(x))
# define omp_get_max_threads()   1
# define omp_get_thread_num()    0
# define omp_get_num_threads()   1
# define omp_lock_t              size_t
# define omp_init_lock(x)        ((void)(x))
# define omp_destroy_lock(x)     ((void)(x))
# define omp_set_lock(x)         ((void)(x))
# define omp_unset_lock(x)       ((void)(x))
#endif

#endif /* __OMP_WRAPPER_H__ */
