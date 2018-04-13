############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(moleculo CXX)


add_executable(spades-truseq-scfcorrection
               main.cpp)
target_link_libraries(spades-truseq-scfcorrection common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-truseq-scfcorrection PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-truseq-scfcorrection
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/scaffold_correction"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")
