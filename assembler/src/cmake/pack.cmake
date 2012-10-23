# -*- cmake -*-

# Binary stuff
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
if (APPLE)
  set(CPACK_GENERATOR "STGZ;TGZ;TBZ2;PackageMaker")
else()
  set(CPACK_GENERATOR "STGZ;TGZ;TBZ2")
endif()

set(CPACK_PACKAGE_NAME "SPAdes")
set(CPACK_PACKAGE_VENDOR "Saint Petersburg Academic University")
set(CPACK_PACKAGE_VERSION "2.3.0")
set(CPACK_PACKAGE_VERSION_MAJOR "2")
set(CPACK_PACKAGE_VERSION_MINOR "3")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_STRIP_FILES bin/debruijn bin/hammer)

# Source stuff
set(CPACK_SOURCE_GENERATOR "TBZ2")
set(CPACK_SOURCE_IGNORE_FILES tools test web_service online_vis cap)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")

include(CPack)
