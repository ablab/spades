cmake_minimum_required (VERSION 2.6)
project (samtools) 

OPTION(BUILD_SHARED_LIBS "Build package with shared libraries." OFF)

add_definitions(-D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE -DMD5SUM_MAIN -D_CURSES_LIB=1)

add_library (bam
# bgzf.c
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
 knetfile.c
 bam_sort.c
 sam_header.c
 bam_reheader.c
 kprobaln.c
 bam_cat.c) 
#because of bgzf redefinition
target_link_libraries(BamTools)

add_library (samtools
 bam_tview.c
 bam_plcmd.c
 sam_view.c
 bam_rmdup.c
 bam_rmdupse.c
 bam_mate.c
 bam_stat.c
 bam_color.c
# bam_qa.c
 bamtk.c
 padding.c
 kaln.c
 bam2bcf.c
 bam2bcf_indel.c
 errmod.c
 sample.c
 cut_target.c
 phase.c
 bam2depth.c)

target_link_libraries (samtools curses bam m z bcf)


include_directories("./")

add_library(bcf
 bcftools/bcf.c
 bcftools/vcf.c
 bcftools/bcfutils.c
 bcftools/prob1.c
 bcftools/em.c
 bcftools/kfunc.c
 bcftools/kmin.c
 bcftools/index.c
 bcftools/fet.c
 bcftools/mut.c
 bcftools/bcf2qcall.c)


add_executable(bcftools
 bcftools/call1.c
 bcftools/main.c
 kstring.c
 bgzf.c
 knetfile.c
 bedidx.c)

target_link_libraries (bcftools bcf m z)

#install(TARGETS samtools bcftools DESTINATION bin)
#install(TARGETS bam bcf DESTINATION lib)
#install(
# FILES
#  bam.h
#  sam.h
#  khash.h
#  faidx.h
#  bgzf.h
#  sam_header.h
# DESTINATION include)

