//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/reads/io_helper.hpp"

namespace io {

class SequencingLibraryBase;

PairedStream paired_easy_reader(const SequencingLibraryBase &lib,
                                bool followed_by_rc,
                                size_t insert_size,
                                bool use_orientation = true,
                                FileReadFlags flags = FileReadFlags(),
                                ThreadPool::ThreadPool *pool = nullptr);
ReadStreamList<SingleRead> single_easy_readers(const SequencingLibraryBase &lib,
                                               bool followed_by_rc,
                                               bool including_paired_reads,
                                               bool handle_Ns = true,
                                               FileReadFlags flags = FileReadFlags(),
                                               ThreadPool::ThreadPool *pool = nullptr);
SingleStream single_easy_reader(const SequencingLibraryBase &lib,
                                bool followed_by_rc,
                                bool including_paired_reads,
                                bool handle_Ns = true,
                                FileReadFlags flags = FileReadFlags(),
                                ThreadPool::ThreadPool *pool = nullptr);
ReadStreamList<SingleRead> merged_easy_readers(const SequencingLibraryBase &lib,
                                               bool followed_by_rc,
                                               bool handle_Ns = true,
                                               FileReadFlags flags = FileReadFlags(),
                                               ThreadPool::ThreadPool *pool = nullptr);
SingleStream merged_easy_reader(const SequencingLibraryBase &lib,
                                bool followed_by_rc,
                                bool handle_Ns = true,
                                FileReadFlags flags = FileReadFlags(),
                                ThreadPool::ThreadPool *pool = nullptr);

}
