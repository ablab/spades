############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(bwa C)

add_executable(bwa               
               bamlite.c  
               bwa.c
               bwase.c 
               bwtaln.c  
               bwtgap.c  
               bwtindex.c 
               bwt_lite.c
               bwtsw2_aux.c 
               bwtsw2_core.c  
               bwtsw2_pair.c
               fastmap.c
               kstring.c
               main.c   
               simple_dp.c
               utils.c
               bntseq.c
               bwape.c
               bwaseqio.c
               bwt.c  
               bwt_gen.c
               bwtio.c  
               bwtmisc.c
               bwtsw2_chain.c
               bwtsw2_main.c
               cs2nt.c     
               is.c    
               ksw.c   
               QSufSort.c
               stdaln.c)

target_link_libraries(bwa ${COMMON_LIBRARIES} m)

if (SPADES_STATIC_BUILD)
  set_target_properties(bwa PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS bwa
        RUNTIME DESTINATION bin)
