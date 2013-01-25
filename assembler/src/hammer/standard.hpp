//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * standart.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once
#include "standard_base.hpp"

//==omp
#include "openmp_wrapper.h"

//==sys
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/syscall.h>

//==our
// utils
#include "cpp_utils.hpp"
#include "logger/logger.hpp"

// io
#include "io/ireader.hpp"
#include "io/converting_reader_wrapper.hpp"
