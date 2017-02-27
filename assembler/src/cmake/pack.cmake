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
set(CPACK_PACKAGE_DESCRIPTION_FILE "${SPADES_MAIN_SRC_DIR}/../README")
set(CPACK_RESOURCE_FILE_LICENSE "${SPADES_MAIN_SRC_DIR}/../LICENSE")
set(CPACK_PACKAGE_VERSION "3.10.1")
set(CPACK_PACKAGE_VERSION_MAJOR "3")
set(CPACK_PACKAGE_VERSION_MINOR "10")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_STRIP_FILES bin/spades bin/hammer bin/ionhammer bin/dipspades bin/spades-bwa bin/corrector bin/scaffold_correction)

# Source stuff
set(CPACK_SOURCE_GENERATOR "TBZ2")
set(CPACK_SOURCE_IGNORE_FILES tools test web_service online_vis cap)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")

include(CPack)
