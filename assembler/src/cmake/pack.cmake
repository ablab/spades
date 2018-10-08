# -*- cmake -*-

# Binary stuff
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
if (APPLE)
  set(CPACK_GENERATOR "STGZ;TGZ;PackageMaker")
else()
  set(CPACK_GENERATOR "STGZ;TGZ")
endif()

set(CPACK_PACKAGE_NAME "SPAdes")
set(CPACK_PACKAGE_VENDOR "Saint Petersburg State University")
#set(CPACK_PACKAGE_DESCRIPTION_FILE "${SPADES_MAIN_SRC_DIR}/../README")
set(CPACK_RESOURCE_FILE_LICENSE "${SPADES_MAIN_SRC_DIR}/../LICENSE")
set(CPACK_PACKAGE_VERSION "3.13.0")
set(CPACK_PACKAGE_VERSION_MAJOR "3")
set(CPACK_PACKAGE_VERSION_MINOR "13")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_STRIP_FILES bin/spades-bwa bin/spades-core bin/spades-corrector-core bin/spades-gbuilder bin/spades-gmapper bin/spades-hammer bin/spades-ionhammer bin/spades-kmercount bin/spades-truseq-scfcorrection)


# Source stuff
set(CPACK_SOURCE_GENERATOR "TBZ2")
set(CPACK_SOURCE_IGNORE_FILES tools test web_service online_vis cap)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")

include(CPack)
