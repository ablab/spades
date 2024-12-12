//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/**
 * @file    parser.cpp
 * @author  Mariya Fomkina
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * Parser is the parent class for all streams that read data from
 * different file types (fastq, fasta, sam etc).
 * This file contains functions that are used to select exact parser
 * according to extension.
 */

#include "config.hpp"
#include "parser.hpp"

#include "file_read_flags.hpp"
#include "fasta_fastq_gz_parser.hpp"
#include "fasta_fastq_zstd_parser.hpp"
#include "io/sam/bam_parser.hpp"
#ifdef SPADES_USE_NCBISDK
# include "io/sra/sra_parser.hpp"
#endif

namespace io {

/*
 * Select parser type according to file extension.
 *
 * @param filename The name of the file to be opened.
 * @param offset The offset of the read quality.

 * @return Pointer to the new parser object with these filename and
 * offset.
 */
Parser* SelectParser(const std::filesystem::path& filename,
                     FileReadFlags flags) {
  if (filename.extension() == ".bam")
      return new BAMParser(filename, flags);
#ifdef SPADES_USE_NCBISDK
  else if (filename.extension() == ".sra")
      return new SRAParser(filename, flags);
#endif
#ifdef SPADES_USE_ZSTD
  else if (filename.extension() == ".zst")
      return new FastaFastqZstdParser(filename, flags);
#endif

  return new FastaFastqGzParser(filename, flags);
}

}
