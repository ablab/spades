# -*- cmake -*-

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(samtools C)

add_definitions(-D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE -DMD5SUM_MAIN -D_CURSES_LIB=1)

add_library (samtools
  kstring.c
  bam_aux.c
  bam.c
  bam_import.c
  sam.c
  bam_index.c
  bam_pileup.c
  bam_lpileup.c
  bam_md.c
  razf.c
  faidx.c
  bedidx.c
  bam_sort.c
  sam_header.c
  bam_reheader.c
  kprobaln.c
  bam_cat.c)
target_link_libraries(samtools BamTools)
