#==================================================================================================
#
#   Alternative CMake file implementing function for importing Blaze in CMake based projects.
#
#   Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
#
#   This file is part of the Blaze library. You can redistribute it and/or modify it under
#   the terms of the New (Revised) BSD License. Redistribution and use in source and binary
#   forms, with or without modification, are permitted provided that the following conditions
#   are met:
#
#   1. Redistributions of source code must retain the above copyright notice, this list of
#      conditions and the following disclaimer.
#   2. Redistributions in binary form must reproduce the above copyright notice, this list
#      of conditions and the following disclaimer in the documentation and/or other materials
#      provided with the distribution.
#   3. Neither the names of the Blaze development group nor the names of its contributors
#      may be used to endorse or promote products derived from this software without specific
#      prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
#   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
#   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
#   SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
#   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
#   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
#   DAMAGE.
#
#==================================================================================================

set(Blaze_Import_FILE_PATH "${CMAKE_CURRENT_LIST_FILE}")
set(Blaze_Import_DIR_PATH "${CMAKE_CURRENT_LIST_DIR}")

function(Blaze_Import)

   macro(msg m)
      if(Blaze_Import_DEBUG OR NOT Blaze_Import_QUIET)
         message(STATUS ${m})
      endif()
   endmacro()

   macro(msg_db m)
      if(Blaze_Import_DEBUG)
         message(STATUS "DEBUG: ${m}")
      endif()
   endmacro()

   cmake_minimum_required(VERSION 3.5)

   #==================================================================================================
   #   Input Parsing
   #==================================================================================================

      set(options
         REQUIRED  # Throw error if Blaze can't be imported
         HELP      # Print Help message
         QUIET     # Suppress messages if not in Debug mode
         DEBUG     # Print additional information messages
      )
      set(oneValueArgs
         BLAS                    # ON/OFF
         BLAS_64BIT              # ON/OFF
         BLAS_PARALLEL           # ON/OFF
         BLAS_MV                 # ON/OFF
         BLAS_MM                 # ON/OFF
         BLAS_INCLUDE            # <cblas.h>/<mkl_blas.h> or another blas header
         LAPACK                  # ON/OFF
         THREADING               # HPX/C++11/Boost/OpenMP/off
         CACHE_SIZE              # auto or an unsinged long (UL) value
         VECTORIZATION           # ON/OFF
         TRANSPOSE_FLAG          # columnVector/rowVector
         STORAGE_ORDER           # rowMajor/ColumnMajor
         PADDING                 # ON/OFF
         STREAMING               # ON/OFF
         OPTIMIZED_KERNELS       # ON/OFF
         DEFAULT_INITIALIZATION  # ON/OFF
         STRONG_INLINE           # ON/OFF
         ALWAYS_INLINE           # ON/OFF
         RESTRICT                # ON/OFF
         DEBUG_MODE              # ON/OFF
         FUNCTION_TRACES         # ON/OFF
         INTERNAL_ASSERTION      # ON/OFF
         USER_ASSERTION          # ON/OFF
         MPI_PARALLEL_MODE       # ON/OFF
         RANDOM                  # (Pseudo) Random Number Generator, like std::mt19937
         # All THRESHOLDS are unsinged long (UL) value
         THRESHOLD_DMATDVECMULT
         THRESHOLD_TDMATDVECMULT
         THRESHOLD_TDVECDMATMULT
         THRESHOLD_TDVECTDMATMULT
         THRESHOLD_DMATDMATMULT
         THRESHOLD_DMATTDMATMULT
         THRESHOLD_TDMATDMATMULT
         THRESHOLD_TDMATTDMATMULT
         THRESHOLD_DMATSMATMULT
         THRESHOLD_TDMATSMATMULT
         THRESHOLD_TSMATDMATMULT
         THRESHOLD_TSMATTDMATMULT
         THRESHOLD_SMP_DVECASSIGN
         THRESHOLD_SMP_DVECSCALARMULT
         THRESHOLD_SMP_DVECDVECADD
         THRESHOLD_SMP_DVECDVECSUB
         THRESHOLD_SMP_DVECDVECMULT
         THRESHOLD_SMP_DVECDVECDIV
         THRESHOLD_SMP_DVECDVECOUTER
         THRESHOLD_SMP_DMATDVECMULT
         THRESHOLD_SMP_TDMATDVECMULT
         THRESHOLD_SMP_TDVECDMATMULT
         THRESHOLD_SMP_TDVECTDMATMULT
         THRESHOLD_SMP_DMATSVECMULT
         THRESHOLD_SMP_TDMATSVECMULT
         THRESHOLD_SMP_TSVECDMATMULT
         THRESHOLD_SMP_TSVECTDMATMULT
         THRESHOLD_SMP_SMATDVECMULT
         THRESHOLD_SMP_TSMATDVECMULT
         THRESHOLD_SMP_TDVECSMATMULT
         THRESHOLD_SMP_TDVECTSMATMULT
         THRESHOLD_SMP_SMATSVECMULT
         THRESHOLD_SMP_TSMATSVECMULT
         THRESHOLD_SMP_TSVECSMATMULT
         THRESHOLD_SMP_TSVECTSMATMULT
         THRESHOLD_SMP_DMATASSIGN
         THRESHOLD_SMP_DMATSCALARMULT
         THRESHOLD_SMP_DMATDMATADD
         THRESHOLD_SMP_DMATTDMATADD
         THRESHOLD_SMP_DMATDMATSUB
         THRESHOLD_SMP_DMATTDMATSUB
         THRESHOLD_SMP_DMATDMATSCHUR
         THRESHOLD_SMP_DMATTDMATSCHUR
         THRESHOLD_SMP_DMATDMATMULT
         THRESHOLD_SMP_DMATTDMATMULT
         THRESHOLD_SMP_TDMATDMATMULT
         THRESHOLD_SMP_TDMATTDMATMULT
         THRESHOLD_SMP_DMATSMATMULT
         THRESHOLD_SMP_DMATTSMATMULT
         THRESHOLD_SMP_TDMATSMATMULT
         THRESHOLD_SMP_TDMATTSMATMULT
         THRESHOLD_SMP_SMATDMATMULT
         THRESHOLD_SMP_SMATTDMATMULT
         THRESHOLD_SMP_TSMATDMATMULT
         THRESHOLD_SMP_TSMATTDMATMULT
         THRESHOLD_SMP_SMATSMATMULT
         THRESHOLD_SMP_SMATTSMATMULT
         THRESHOLD_SMP_TSMATSMATMULT
         THRESHOLD_SMP_TSMATTSMATMULT
         THRESHOLD_SMP_DMATREDUCE
         THRESHOLD_SMP_SMATREDUCE
      )
      set(multiValueArgs )

      cmake_parse_arguments(Blaze_Import
         "${options}"
         "${oneValueArgs}"
         "${multiValueArgs}"
         ${ARGN}
      )

   #==================================================================================================
   #   Include Guard
   #==================================================================================================

      if(blaze_FOUND)
         if(NOT Blaze_Import_QUIET OR Blaze_Import_DEBUG)
            message(WARNING "Blaze was already imported and will not be imported again")
         endif()
      else()

   #==================================================================================================
   #   Blaze Version
   #==================================================================================================

      file(READ ${Blaze_Import_DIR_PATH}/../blaze/system/Version.h BLAZE_VERSION_FILE)

      string(REGEX MATCH "#define BLAZE_MAJOR_VERSION ([0-9]*)" _BLAZE_MAJOR_VERSION ${BLAZE_VERSION_FILE})
      set(BLAZE_MAJOR_VERSION ${CMAKE_MATCH_1})

      string(REGEX MATCH "#define BLAZE_MINOR_VERSION ([0-9]*)" _BLAZE_MINOR_VERSION ${BLAZE_VERSION_FILE})
      set(BLAZE_MINOR_VERSION ${CMAKE_MATCH_1})

      set(BLAZE_VERSION ${BLAZE_MAJOR_VERSION}.${BLAZE_MINOR_VERSION})

      msg("Configuring blaze version ${BLAZE_VERSION}")

   #==================================================================================================
   #   Help Message
   #==================================================================================================

      if(Blaze_Import_HELP)
         message("\n\t\t______________________________")
         message("\t\t\tBlaze_Import Module Info")
         message("\t\t______________________________\n")

         message("Description")

         message("________________________________________")

      endif()

   #==================================================================================================
   #   Library Definition
   #==================================================================================================

      add_library(Blaze INTERFACE)
      target_include_directories(Blaze SYSTEM INTERFACE
         $<BUILD_INTERFACE:${Blaze_Import_DIR_PATH}/../>
         $<INSTALL_INTERFACE:include>
         )

      if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
         message(STATUS "Blaze requires compiling in mode supporting C++14 features.")
      else()
         target_compile_features(Blaze INTERFACE cxx_std_14)
      endif()

   #==================================================================================================
   #   Shared Memory Parallelization
   #==================================================================================================

      if(Blaze_Import_THREADING  STREQUAL "HPX")
         find_package(HPX REQUIRED)
         target_compile_definitions(blaze INTERFACE BLAZE_USE_HPX_THREADS)
         target_include_directories(blaze INTERFACE ${HPX_INCLUDE_DIRS})
         target_link_libraries(blaze INTERFACE ${HPX_LIBRARIES})
         msg("Configured HPX for multithreading.")
      elseif(Blaze_Import_THREADING  STREQUAL "C++11")
         find_package(Threads REQUIRED)
         target_compile_definitions(Blaze INTERFACE BLAZE_USE_CPP_THREADS)
         target_link_libraries(Blaze INTERFACE ${CMAKE_THREAD_LIBS_INIT})
         msg("Configured C++11 for multithreading.")
      elseif(Blaze_Import_THREADING  STREQUAL "BOOST")
         find_package(Boost REQUIRED COMPONENTS thread)
         target_compile_definitions(blaze INTERFACE BLAZE_USE_BOOST_THREADS)
         target_include_directories(blaze INTERFACE ${Boost_INCLUDE_DIRS})
         target_compile_definitions(blaze INTERFACE BOOST_ALL_DYN_LINK ) # Needed for Visual Studio 2015
         target_link_libraries(blaze INTERFACE ${Boost_LIBRARIES})
         msg("Configured BOOST for multithreading.")
      elseif(Blaze_Import_THREADING  STREQUAL "OpenMP")
         find_package(OpenMP)
         if (OPENMP_FOUND)
            target_compile_options(blaze INTERFACE ${OpenMP_CXX_FLAGS})
            if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
               target_link_libraries(blaze INTERFACE ${OpenMP_CXX_FLAGS}) # Needed for GNU and Clang
            endif ()
         msg("Configured OpenMP for multithreading.")
         else ()
            message(WARNING "OpenMP not found. Blaze is running on a single thread. Try C++11 or Boost.")
         endif ()
      elseif(Blaze_Import_THREADING  STREQUAL "off")
         msg("Multithreading for Blaze is turned off.")
      elseif("${Blaze_Import_THREADING}"  STREQUAL "")
         msg_db("Multithreading not configured. Using Blaze's defaults.")
      endif()

   #==================================================================================================
   #   BLAS
   #==================================================================================================

      if(Blaze_Import_BLAS)
         target_compile_definitions( Blaze INTERFACE BLAZE_BLAS_MODE=1 )
         find_package(BLAS REQUIRED)
         target_link_libraries(Blaze INTERFACE $<BUILD_INTERFACE:${BLAS_LIBRARIES}>)
         target_compile_options(Blaze INTERFACE $<BUILD_INTERFACE:${BLAS_LINKER_FLAGS}>)
         msg("Configuring BLAS : ON")
      elseif("${Blaze_Import_BLAS}" STREQUAL "")
         msg_db("Using default configuration for BLAS.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_BLAS_MODE=0 )
         msg("Configuring BLAS : OFF")
      endif()

      if(Blaze_Import_BLAS_64BIT)
         target_compile_definitions( Blaze INTERFACE BLAZE_BLAS_IS_64BIT=1 )
         msg_db("Configuring 64-bit BLAS : ON")
      elseif("${Blaze_Import_BLAS_64BIT}" STREQUAL "")
         msg_db("Using default configuration for 64-bit BLAS.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_BLAS_IS_64BIT=0 )
         msg_db("Configuring 64-bit BLAS : OFF")
      endif()

      if(Blaze_Import_BLAS_PARALLEL)
         target_compile_definitions( Blaze INTERFACE BLAZE_BLAS_IS_PARALLEL=1 )
         msg_db("Configuring parallel BLAS : ON")
      elseif("${Blaze_Import_BLAS_PARALLEL}" STREQUAL "")
         msg_db("Using default configuration for parallel BLAS.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_BLAS_IS_PARALLEL=0 )
         msg_db("Configuring parallel BLAS : OFF")
      endif()

      if(Blaze_Import_BLAS_MV)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION=1 )
         msg_db("Configuring BLAS Matrix/Vector Multiplication : ON")
      elseif("${Blaze_Import_BLAS_MV}" STREQUAL "")
         msg_db("Using default configuration for BLAS Matrix/Vector Multiplication.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION=0 )
         msg_db("Configuring BLAS Matrix/Vector Multiplication : OFF")
      endif()

      if(Blaze_Import_BLAS_MM)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION=1 )
         msg_db("Configuring BLAS Matrix/Matrix Multiplication : ON")
      elseif("${Blaze_Import_BLAS_MM}" STREQUAL "")
         msg_db("Using default configuration for BLAS Matrix/Matrix Multiplication.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION=0 )
         msg_db("Configuring BLAS Matrix/Matrix Multiplication : OFF")
      endif()

      if(Blaze_Import_BLAS_INCLUDE)
         target_compile_definitions( Blaze INTERFACE BLAZE_BLAS_INCLUDE_FILE=${Blaze_Import_BLAS_INCLUDE} )
         msg_db("Configuring BLAS Include File : ON")
      elseif("${Blaze_Import_BLAS_INCLUDE}" STREQUAL "")
         msg_db("Using default configuration for BLAS Include File.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_BLAS_INCLUDE_FILE=<cblas.h> )
         msg_db("Configuring BLAS Include File : OFF")
      endif()

   #==================================================================================================
   #   LAPACK
   #==================================================================================================

      if(${Blaze_Import_LAPACK})
         find_package(LAPACK REQUIRED)
         target_link_libraries(Blaze INTERFACE $<BUILD_INTERFACE:${LAPACK_LIBRARIES}>)
         target_compile_options(Blaze INTERFACE $<BUILD_INTERFACE:${LAPACK_LINKER_FLAGS}>)
         msg("Configuring LAPACK : ON")
      elseif("${Blaze_Import_LAPACK}" STREQUAL "")
         msg_db("Using default configuration for LAPACK.")
      else()
         msg("Configuring LAPACK : OFF")
      endif()

   #==================================================================================================
   #   Vectorization
   #==================================================================================================

      if(Blaze_Import_VECTORIZATION)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_VECTORIZATION=1 )
         msg("Configuring Vectorization : ON")
      elseif("${Blaze_Import_VECTORIZATION}" STREQUAL "")
         msg_db("Using default configuration for Vectorization.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_VECTORIZATION=0 )
         msg("Configuring Vectorization : OFF")
      endif()

   #==================================================================================================
   #   Cache
   #==================================================================================================

      if ("${Blaze_Import_CACHE_SIZE}" STREQUAL "auto")
         msg("Automatic Cache Size Configuration")
         set(flag 1)
         if (WIN32)
           execute_process(COMMAND wmic cpu get L3CacheSize
                       OUTPUT_VARIABLE tmp
                       RESULT_VARIABLE flag
                       ERROR_QUIET)
           if (flag)
             execute_process(COMMAND wmic cpu get L2CacheSize
                         OUTPUT_VARIABLE tmp
                         RESULT_VARIABLE flag
                         ERROR_QUIET)
           endif (flag)
           if (flag)
             execute_process(COMMAND wmic cpu get L1CacheSize
                         OUTPUT_VARIABLE tmp
                         RESULT_VARIABLE flag
                         ERROR_QUIET)
           endif (flag)
         endif (WIN32)

         if (UNIX)
           execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index3/size
                       OUTPUT_VARIABLE tmp
                       RESULT_VARIABLE flag
                       ERROR_QUIET)
           if (flag)
             execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index2/size
                         OUTPUT_VARIABLE tmp
                         RESULT_VARIABLE flag
                         ERROR_QUIET)
           endif (flag)
           if (flag)
             execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index1/size
                         OUTPUT_VARIABLE tmp
                         RESULT_VARIABLE flag
                         ERROR_QUIET)
           endif (flag)
         endif (UNIX)

         if (APPLE)
           execute_process(COMMAND sysctl -n hw.l3cachesize
                       OUTPUT_VARIABLE tmp
                       RESULT_VARIABLE flag
                       ERROR_QUIET)
           if (flag OR NOT tmp)
             execute_process(COMMAND sysctl -n hw.l2cachesize
                         OUTPUT_VARIABLE tmp
                         RESULT_VARIABLE flag
                         ERROR_QUIET)
           endif ()
           if (flag OR NOT tmp)
             execute_process(COMMAND sysctl -n hw.l1icachesize
                         OUTPUT_VARIABLE tmp
                         RESULT_VARIABLE flag
                         ERROR_QUIET)
           endif ()

           if (flag EQUAL 0)
             math(EXPR tmp ${tmp}/1024)  # If successful convert to kibibytes to comply with rest
           endif (flag EQUAL 0)
         endif (APPLE)

         if (flag)
           message(WARNING "Cache size not found automatically. Using default value as cache size.")
           set(tmp "3072")
         endif (flag)

         string(REGEX MATCH "([0-9][0-9]+)" tmp ${tmp}) # Get a number containing at least 2 digits in the string tmp
         math(EXPR BLAZE_CACHE_SIZE ${tmp}*1024) # Convert to bytes (assuming that the value is given in kibibytes)
         set(BLAZE_CACHE_SIZE "${BLAZE_CACHE_SIZE}UL")
         target_compile_definitions( Blaze INTERFACE BLAZE_CACHE_SIZE=${BLAZE_CACHE_SIZE} )
         msg("Configuring Cache Size : ${BLAZE_CACHE_SIZE}")
      elseif(${Blaze_Import_CACHE_SIZE})
         target_compile_definitions( Blaze INTERFACE BLAZE_CACHE_SIZE=${Blaze_Import_CACHE_SIZE} )
         msg("Configuring Cache Size : ${Blaze_Import_CACHE_SIZE}")
      else()
         msg_db("Using default Cache Size.")
      endif ()

   #==================================================================================================
   #   Transpose Flag
   #==================================================================================================

      if(Blaze_Import_TRANSPOSE_FLAG STREQUAL "columnVector")
         target_compile_definitions( Blaze INTERFACE BLAZE_DEFAULT_TRANSPOSE_FLAG=blaze::columnVector )
         msg("Configuring Transpose Flag : columnVector")
      elseif(Blaze_Import_TRANSPOSE_FLAG STREQUAL "rowVector")
         target_compile_definitions( Blaze INTERFACE BLAZE_DEFAULT_TRANSPOSE_FLAG=blaze::rowVector )
         msg("Configuring Transpose Flag : rowVector")
      else("${Blaze_Import_TRANSPOSE_FLAG}" STREQUAL "")
         msg_db("Using default configuration for Transpose Flag.")
      endif()

   #==================================================================================================
   #   Storage Order
   #==================================================================================================

      if(Blaze_Import_STORAGE_ORDER STREQUAL "rowMajor")
         target_compile_definitions( Blaze INTERFACE BLAZE_DEFAULT_STORAGE_ORDER=blaze::rowMajor )
         msg("Configuring Storage Order : rowMajor")
      elseif(Blaze_Import_STORAGE_ORDER STREQUAL "columnMajor")
         target_compile_definitions( Blaze INTERFACE BLAZE_DEFAULT_STORAGE_ORDER=blaze::columnMajor )
         msg("Configuring Storage Order : columnMajor")
      else("${Blaze_Import_STORAGE_ORDER}" STREQUAL "")
         msg_db("Using default configuration for Storage Order.")
      endif()

   #==================================================================================================
   #   Optimization
   #==================================================================================================

      if(Blaze_Import_PADDING)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_PADDING=1 )
         msg("Configuring Padding : ON")
      elseif("${Blaze_Import_PADDING}" STREQUAL "")
         msg_db("Using default configuration for Padding.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_PADDING=0 )
         msg("Configuring Padding : OFF")
      endif()

      if(Blaze_Import_STREAMING)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_STREAMING=1 )
         msg("Configuring Streaming : ON")
      elseif("${Blaze_Import_STREAMING}" STREQUAL "")
         msg_db("Using default configuration for Streaming.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_STREAMING=0 )
         msg("Configuring Streaming : OFF")
      endif()

      if(Blaze_Import_OPTIMIZED_KERNELS)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_OPTIMIZED_KERNELS=1 )
         msg("Configuring Optimized Kernels : ON")
      elseif("${Blaze_Import_OPTIMIZED_KERNELS}" STREQUAL "")
         msg_db("Using default configuration for Optimized Kernels.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_OPTIMIZED_KERNELS=0 )
         msg("Configuring Optimized Kernels : OFF")
      endif()

      if(Blaze_Import_DEFAULT_INITIALIZATION)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_DEFAULT_INITIALIZATION=1 )
         msg("Configuring Default Initialization : ON")
      elseif("${Blaze_Import_DEFAULT_INITIALIZATION}" STREQUAL "")
         msg_db("Using default configuration for Default Initialization.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_DEFAULT_INITIALIZATION=0 )
         msg("Configuring Default Initialization : OFF")
      endif()

   #==================================================================================================
   #   Inlining
   #==================================================================================================

      if(Blaze_Import_DEFAULT_STRONG_INLINE)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_STRONG_INLINE=1 )
         msg("Configuring Strong Inlining : ON")
      elseif("${Blaze_Import_DEFAULT_STRONG_INLINE}" STREQUAL "")
         msg_db("Using default configuration for Strong Inlining.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_STRONG_INLINE=0 )
         msg("Configuring Strong Inlining : OFF")
      endif()

      if(Blaze_Import_DEFAULT_ALWAYS_INLINE)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_ALWAYS_INLINE=1 )
         msg("Configuring Always Inlining : ON")
      elseif("${Blaze_Import_DEFAULT_ALWAYS_INLINE}" STREQUAL "")
         msg_db("Using default configuration for Always Inlining.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_ALWAYS_INLINE=0 )
         msg("Configuring Always Inlining : OFF")
      endif()

   #==================================================================================================
   #   Thresholds
   #==================================================================================================

      if(Blaze_Import_THRESHOLD_DMATDVECMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_DMATDVECMULT_THRESHOLD=${Blaze_Import_THRESHOLD_DMATDVECMULT} )
         msg_db("Configuring Row-Major Dense Matrix/Vector Multiplication Threshold : ${Blaze_Import_THRESHOLD_DMATDVECMULT}")
      else()
         msg_db("Using default configuration for Row-Major Dense Matrix/Vector Multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_TDMATDVECMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_TDMATDVECMULT_THRESHOLD=${Blaze_Import_THRESHOLD_TDMATDVECMULT} )
         msg_db("Configuring Column-Major Dense Matrix/Vector Multiplication Threshold : ${Blaze_Import_THRESHOLD_TDMATDVECMULT}")
      else()
         msg_db("Using default configuration for Column-Major Dense Matrix/Vector Multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_TDVECDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_TDVECDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_TDVECDMATMULT} )
         msg_db("Configuring Dense Vector/row-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_TDVECDMATMULT}")
      else()
         msg_db("Using default configuration for Dense Vector/row-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_TDVECTDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_TDVECTDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_TDVECTDMATMULT} )
         msg_db("Configuring Dense Vector/column-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_TDVECTDMATMULT}")
      else()
         msg_db("Using default configuration for Dense Vector/column-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_DMATDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_DMATDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_DMATDMATMULT} )
         msg_db("Configuring Row-major dense matrix/row-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_DMATDMATMULT}")
      else()
         msg_db("Using default configuration for Row-major dense matrix/row-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_DMATTDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_DMATTDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_DMATTDMATMULT} )
         msg_db("Configuring Row-major dense matrix/column-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_DMATTDMATMULT}")
      else()
         msg_db("Using default configuration for Row-major dense matrix/column-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_TDMATDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_TDMATDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_TDMATDMATMULT} )
         msg_db("Configuring Column-major dense matrix/row-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_TDMATDMATMULT}")
      else()
         msg_db("Using default configuration for Column-major dense matrix/row-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_TDMATTDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_TDMATTDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_TDMATTDMATMULT} )
         msg_db("Configuring Column-major dense matrix/column-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_TDMATTDMATMULT}")
      else()
         msg_db("Using default configuration for Column-major dense matrix/column-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_DMATSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_DMATSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_DMATSMATMULT} )
         msg_db("Configuring Row-major dense matrix/row-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_DMATSMATMULT}")
      else()
         msg_db("Using default configuration for Row-major dense matrix/row-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_TDMATSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_TDMATSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_TDMATSMATMULT} )
         msg_db("Configuring Column-major dense matrix/row-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_TDMATSMATMULT}")
      else()
         msg_db("Using default configuration for Column-major dense matrix/row-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_TSMATDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_TSMATDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_TSMATDMATMULT} )
         msg_db("Configuring Column-major sparse matrix/row-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_TSMATDMATMULT}")
      else()
         msg_db("Using default configuration for Column-major sparse matrix/row-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_TSMATTDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_TSMATTDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_TSMATTDMATMULT} )
         msg_db("Configuring Column-major sparse matrix/column-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_TSMATTDMATMULT}")
      else()
         msg_db("Using default configuration for Column-major sparse matrix/column-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DVECASSIGN)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DVECASSIGN_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DVECASSIGN} )
         msg_db("Configuring SMP dense vector assignment Threshold : ${Blaze_Import_THRESHOLD_SMP_DVECASSIGN}")
      else()
         msg_db("Using default configuration for SMP dense vector assignment Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DVECSCALARMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DVECSCALARMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DVECSCALARMULT} )
         msg_db("Configuring SMP dense vector/scalar multiplication/division Threshold : ${Blaze_Import_THRESHOLD_SMP_DVECSCALARMULT}")
      else()
         msg_db("Using default configuration for SMP dense vector/scalar multiplication/division Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DVECDVECADD)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DVECDVECADD_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DVECDVECADD} )
         msg_db("Configuring SMP dense vector/dense vector addition Threshold : ${Blaze_Import_THRESHOLD_SMP_DVECDVECADD}")
      else()
         msg_db("Using default configuration for SMP dense vector/dense vector addition Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DVECDVECSUB)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DVECDVECSUB_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DVECDVECSUB} )
         msg_db("Configuring SMP dense vector/dense vector subtraction Threshold : ${Blaze_Import_THRESHOLD_SMP_DVECDVECSUB}")
      else()
         msg_db("Using default configuration for SMP dense vector/dense vector subtraction Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DVECDVECMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DVECDVECMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DVECDVECMULT} )
         msg_db("Configuring SMP dense vector/dense vector multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_DVECDVECMULT}")
      else()
         msg_db("Using default configuration for SMP dense vector/dense vector multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DVECDVECDIV)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DVECDVECDIV_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DVECDVECDIV} )
         msg_db("Configuring SMP dense vector/dense vector division Threshold : ${Blaze_Import_THRESHOLD_SMP_DVECDVECDIV}")
      else()
         msg_db("Using default configuration for SMP dense vector/dense vector division Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DVECDVECOUTER)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DVECDVECOUTER_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DVECDVECOUTER} )
         msg_db("Configuring SMP dense vector/dense vector outer product Threshold : ${Blaze_Import_THRESHOLD_SMP_DVECDVECOUTER}")
      else()
         msg_db("Using default configuration for SMP dense vector/dense vector outer product Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATDVECMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATDVECMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATDVECMULT} )
         msg_db("Configuring SMP row-major dense matrix/dense vector multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATDVECMULT}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/dense vector multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TDMATDVECMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TDMATDVECMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TDMATDVECMULT} )
         msg_db("Configuring SMP column-major dense matrix/dense vector multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TDMATDVECMULT}")
      else()
         msg_db("Using default configuration for SMP column-major dense matrix/dense vector multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TDVECDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TDVECDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TDVECDMATMULT} )
         msg_db("Configuring SMP dense vector/row-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TDVECDMATMULT}")
      else()
         msg_db("Using default configuration for SMP dense vector/row-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TDVECTDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TDVECTDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TDVECTDMATMULT} )
         msg_db("Configuring SMP dense vector/column-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TDVECTDMATMULT}")
      else()
         msg_db("Using default configuration for SMP dense vector/column-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATSVECMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATSVECMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATSVECMULT} )
         msg_db("Configuring SMP row-major dense matrix/sparse vector multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATSVECMULT}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/sparse vector multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TDMATSVECMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TDMATSVECMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TDMATSVECMULT} )
         msg_db("Configuring SMP column-major dense matrix/sparse vector multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TDMATSVECMULT}")
      else()
         msg_db("Using default configuration for SMP column-major dense matrix/sparse vector multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TSVECDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TSVECDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TSVECDMATMULT} )
         msg_db("Configuring SMP sparse vector/row-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TSVECDMATMULT}")
      else()
         msg_db("Using default configuration for SMP sparse vector/row-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TSVECTDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TSVECTDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TSVECTDMATMULT} )
         msg_db("Configuring SMP sparse vector/column-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TSVECTDMATMULT}")
      else()
         msg_db("Using default configuration for SMP sparse vector/column-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_SMATDVECMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_SMATDVECMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_SMATDVECMULT} )
         msg_db("Configuring SMP row-major sparse matrix/dense vector multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_SMATDVECMULT}")
      else()
         msg_db("Using default configuration for SMP row-major sparse matrix/dense vector multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TSMATDVECMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TSMATDVECMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TSMATDVECMULT} )
         msg_db("Configuring SMP column-major sparse matrix/dense vector multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TSMATDVECMULT}")
      else()
         msg_db("Using default configuration for SMP column-major sparse matrix/dense vector multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TDVECSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TDVECSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TDVECSMATMULT} )
         msg_db("Configuring SMP dense vector/row-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TDVECSMATMULT}")
      else()
         msg_db("Using default configuration for SMP dense vector/row-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TDVECTSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TDVECTSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TDVECTSMATMULT} )
         msg_db("Configuring SMP dense vector/column-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TDVECTSMATMULT}")
      else()
         msg_db("Using default configuration for SMP dense vector/column-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_SMATSVECMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_SMATSVECMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_SMATSVECMULT} )
         msg_db("Configuring SMP row-major sparse matrix/sparse vector multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_SMATSVECMULT}")
      else()
         msg_db("Using default configuration for SMP row-major sparse matrix/sparse vector multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TSMATSVECMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TSMATSVECMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TSMATSVECMULT} )
         msg_db("Configuring SMP column-major sparse matrix/sparse vector multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TSMATSVECMULT}")
      else()
         msg_db("Using default configuration for SMP column-major sparse matrix/sparse vector multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TSVECSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TSVECSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TSVECSMATMULT} )
         msg_db("Configuring SMP sparse vector/row-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TSVECSMATMULT}")
      else()
         msg_db("Using default configuration for SMP sparse vector/row-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TSVECTSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TSVECTSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TSVECTSMATMULT} )
         msg_db("Configuring SMP sparse vector/column-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TSVECTSMATMULT}")
      else()
         msg_db("Using default configuration for SMP sparse vector/column-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATASSIGN)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATASSIGN_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATASSIGN} )
         msg_db("Configuring SMP dense matrix assignment Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATASSIGN}")
      else()
         msg_db("Using default configuration for SMP dense matrix assignment Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATSCALARMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATSCALARMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATSCALARMULT} )
         msg_db("Configuring SMP dense matrix/scalar multiplication/division Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATSCALARMULT}")
      else()
         msg_db("Using default configuration for SMP dense matrix/scalar multiplication/division Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATDMATADD)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATDMATADD_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATDMATADD} )
         msg_db("Configuring SMP row-major dense matrix/row-major dense matrix addition Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATDMATADD}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/row-major dense matrix addition Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATTDMATADD)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATTDMATADD_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATTDMATADD} )
         msg_db("Configuring SMP row-major dense matrix/column-major dense matrix addition Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATTDMATADD}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/column-major dense matrix addition Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATDMATSUB)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATDMATSUB_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATDMATSUB} )
         msg_db("Configuring SMP row-major dense matrix/row-major dense matrix subtraction Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATDMATSUB}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/row-major dense matrix subtraction Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATTDMATSUB)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATTDMATSUB_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATTDMATSUB} )
         msg_db("Configuring SMP row-major dense matrix/column-major dense matrix subtraction Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATTDMATSUB}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/column-major dense matrix subtraction Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATDMATSCHUR)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATDMATSCHUR_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATDMATSCHUR} )
         msg_db("Configuring SMP row-major dense matrix/row-major dense matrix Schur product Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATDMATSCHUR}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/row-major dense matrix Schur product Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATTDMATSCHUR)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATTDMATSCHUR_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATTDMATSCHUR} )
         msg_db("Configuring SMP row-major dense matrix/column-major dense matrix Schur product Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATTDMATSCHUR}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/column-major dense matrix Schur product Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATDMATMULT} )
         msg_db("Configuring SMP row-major dense matrix/row-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATDMATMULT}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/row-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATTDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATTDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATTDMATMULT} )
         msg_db("Configuring SMP row-major dense matrix/column-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATTDMATMULT}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/column-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TDMATDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TDMATDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TDMATDMATMULT} )
         msg_db("Configuring SMP column-major dense matrix/row-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TDMATDMATMULT}")
      else()
         msg_db("Using default configuration for SMP column-major dense matrix/row-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TDMATTDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TDMATTDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TDMATTDMATMULT} )
         msg_db("Configuring SMP column-major dense matrix/column-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TDMATTDMATMULT}")
      else()
         msg_db("Using default configuration for SMP column-major dense matrix/column-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATSMATMULT} )
         msg_db("Configuring SMP row-major dense matrix/row-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATSMATMULT}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/row-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATTSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATTSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATTSMATMULT} )
         msg_db("Configuring SMP row-major dense matrix/column-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATTSMATMULT}")
      else()
         msg_db("Using default configuration for SMP row-major dense matrix/column-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TDMATSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TDMATSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TDMATSMATMULT} )
         msg_db("Configuring SMP column-major dense matrix/row-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TDMATSMATMULT}")
      else()
         msg_db("Using default configuration for SMP column-major dense matrix/row-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TDMATTSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TDMATTSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TDMATTSMATMULT} )
         msg_db("Configuring SMP column-major dense matrix/column-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TDMATTSMATMULT}")
      else()
         msg_db("Using default configuration for SMP column-major dense matrix/column-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_SMATDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_SMATDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_SMATDMATMULT} )
         msg_db("Configuring SMP row-major sparse matrix/row-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_SMATDMATMULT}")
      else()
         msg_db("Using default configuration for SMP row-major sparse matrix/row-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_SMATTDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_SMATTDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_SMATTDMATMULT} )
         msg_db("Configuring SMP row-major sparse matrix/column-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_SMATTDMATMULT}")
      else()
         msg_db("Using default configuration for SMP row-major sparse matrix/column-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TSMATDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TSMATDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TSMATDMATMULT} )
         msg_db("Configuring SMP column-major sparse matrix/row-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TSMATDMATMULT}")
      else()
         msg_db("Using default configuration for SMP column-major sparse matrix/row-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TSMATTDMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TSMATTDMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TSMATTDMATMULT} )
         msg_db("Configuring SMP column-major sparse matrix/column-major dense matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TSMATTDMATMULT}")
      else()
         msg_db("Using default configuration for SMP column-major sparse matrix/column-major dense matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_SMATSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_SMATSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_SMATSMATMULT} )
         msg_db("Configuring SMP row-major sparse matrix/row-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_SMATSMATMULT}")
      else()
         msg_db("Using default configuration for SMP row-major sparse matrix/row-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_SMATTSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_SMATTSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_SMATTSMATMULT} )
         msg_db("Configuring SMP row-major sparse matrix/column-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_SMATTSMATMULT}")
      else()
         msg_db("Using default configuration for SMP row-major sparse matrix/column-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TSMATSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TSMATSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TSMATSMATMULT} )
         msg_db("Configuring SMP column-major sparse matrix/row-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TSMATSMATMULT}")
      else()
         msg_db("Using default configuration for SMP column-major sparse matrix/row-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_TSMATTSMATMULT)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_TSMATTSMATMULT_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_TSMATTSMATMULT} )
         msg_db("Configuring SMP column-major sparse matrix/column-major sparse matrix multiplication Threshold : ${Blaze_Import_THRESHOLD_SMP_TSMATTSMATMULT}")
      else()
         msg_db("Using default configuration for SMP column-major sparse matrix/column-major sparse matrix multiplication Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_DMATREDUCE)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_DMATREDUCE_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_DMATREDUCE} )
         msg_db("Configuring SMP dense matrix reduction Threshold : ${Blaze_Import_THRESHOLD_SMP_DMATREDUCE}")
      else()
         msg_db("Using default configuration for SMP dense matrix reduction Threshold.")
      endif()

      if(Blaze_Import_THRESHOLD_SMP_SMATREDUCE)
         target_compile_definitions( Blaze INTERFACE BLAZE_SMP_SMATREDUCE_THRESHOLD=${Blaze_Import_THRESHOLD_SMP_SMATREDUCE} )
         msg_db("Configuring SMP sparse matrix reduction Threshold : ${Blaze_Import_THRESHOLD_SMP_SMATREDUCE}")
      else()
         msg_db("Using default configuration for SMP sparse matrix reduction Threshold.")
      endif()

   #==================================================================================================
   #   MPI
   #==================================================================================================

      if(Blaze_Import_MPI_PARALLEL_MODE)
         target_compile_definitions( Blaze INTERFACE BLAZE_MPI_PARALLEL_MODE=1 )
         msg("Configuring MPI Parallel Mode : ON")
      elseif("${Blaze_Import_MPI_PARALLEL_MODE}" STREQUAL "")
         msg_db("Using default configuration for MPI Parallel Mode.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_MPI_PARALLEL_MODE=0 )
         msg("Configuring MPI Parallel Mode : OFF")
      endif()

   #==================================================================================================
   #   Random
   #==================================================================================================

      if(Blaze_Import_RANDOM)
         target_compile_definitions( Blaze INTERFACE BLAZE_RANDOM_NUMBER_GENERATOR=${Blaze_Import_RANDOM} )
         msg("Configuring Random Number Generator : ${Blaze_Import_RANDOM}")
      else()
         msg_db("Using default configuration for Random Number Generator.")
      endif()

   #==================================================================================================
   #   Restrict
   #==================================================================================================

      if(Blaze_Import_RESTRICT)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_RESTRICT=1 )
         msg("Configuring Restrict Keyword : ON")
      elseif("${Blaze_Import_RESTRICT}" STREQUAL "")
         msg_db("Using default configuration for Restrict Keyword.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_RESTRICT=0 )
         msg("Configuring Restrict Keyword : OFF")
      endif()

   #==================================================================================================
   #   Debugging
   #==================================================================================================

      if(Blaze_Import_DEBUG_MODE)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_DEBUG_MODE=1 )
         msg("Configuring Debug Mode : ON")
      elseif("${Blaze_Import_DEBUG_MODE}" STREQUAL "")
         if(CMAKE_BUILD_TYPE STREQUAL "Debug")
            target_compile_definitions( Blaze INTERFACE BLAZE_USE_DEBUG_MODE=1 )
            msg_db("Configuring Debug Mode : ON")
         else()
            target_compile_definitions( Blaze INTERFACE BLAZE_USE_DEBUG_MODE=0 )
            msg_db("Configuring Debug Mode : OFF")
         endif()
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_DEBUG_MODE=0 )
         msg("Configuring Debug Mode : OFF")
      endif()

      if(Blaze_Import_FUNCTION_TRACES)
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_FUNCTION_TRACES=1 )
         msg("Configuring Function Traces : ON")
      elseif("${Blaze_Import_FUNCTION_TRACES}" STREQUAL "")
         if(CMAKE_BUILD_TYPE STREQUAL "Debug")
            target_compile_definitions( Blaze INTERFACE BLAZE_USE_FUNCTION_TRACES=1 )
            msg_db("Configuring Function Traces : ON")
         else()
            target_compile_definitions( Blaze INTERFACE BLAZE_USE_FUNCTION_TRACES=0 )
            msg_db("Configuring Function Traces : OFF")
         endif()
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USE_FUNCTION_TRACES=0 )
         msg("Configuring Function Traces : OFF")
      endif()

   #==================================================================================================
   #   Assertions
   #==================================================================================================

      if(Blaze_Import_INTERNAL_ASSERTION)
         target_compile_definitions( Blaze INTERFACE BLAZE_INTERNAL_ASSERTION=1 )
         msg("Configuring Internal Assertion : ON")
      elseif("${Blaze_Import_INTERNAL_ASSERTION}" STREQUAL "")
         msg_db("Using default configuration for Internal Assertion.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_INTERNAL_ASSERTION=0 )
         msg("Configuring Internal Assertion : OFF")
      endif()

      if(Blaze_Import_USER_ASSERTION)
         target_compile_definitions( Blaze INTERFACE BLAZE_USER_ASSERTION=1 )
         msg("Configuring User Assertion : ON")
      elseif("${Blaze_Import_USER_ASSERTION}" STREQUAL "")
         msg_db("Using default configuration for User Assertion.")
      else()
         target_compile_definitions( Blaze INTERFACE BLAZE_USER_ASSERTION=0 )
         msg("Configuring User Assertion : OFF")
      endif()

   #==================================================================================================
   #
   #==================================================================================================

      #set(blaze_FOUND TRUE CACHE BOOL "Indicates whether the blaze library was found successfully.")

      msg("Blaze has been configured successfully")

   endif()


endfunction()
