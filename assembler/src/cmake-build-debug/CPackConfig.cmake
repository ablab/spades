# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


set(CPACK_BUILD_SOURCE_DIRS "/home/yopican/Downloads/bi/git/splitter/assembler/src;/home/yopican/Downloads/bi/git/splitter/assembler/src/cmake-build-debug")
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
set(CPACK_COMPONENTS_ALL "Unspecified;common;runtime;spades;spades-docs;spades-test;spaligner")
set(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
set(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_FILE "/home/yopican/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/233.13135.93/bin/cmake/linux/x64/share/cmake-3.27/Templates/CPack.GenericDescription.txt")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_SUMMARY "Project built using CMake")
set(CPACK_DMG_SLA_USE_RESOURCE_FILE_LICENSE "ON")
set(CPACK_GENERATOR "STGZ;TGZ")
set(CPACK_INNOSETUP_ARCHITECTURE "x64")
set(CPACK_INSTALL_CMAKE_PROJECTS "/home/yopican/Downloads/bi/git/splitter/assembler/src/cmake-build-debug;Project;ALL;/")
set(CPACK_INSTALL_PREFIX "/usr/local")
set(CPACK_MODULE_PATH "/home/yopican/Downloads/bi/git/splitter/assembler/src/cmake;/home/yopican/Downloads/bi/git/splitter/assembler/src/cmake/Modules")
set(CPACK_NSIS_DISPLAY_NAME "SPAdes 3.16.0-dev")
set(CPACK_NSIS_INSTALLER_ICON_CODE "")
set(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
set(CPACK_NSIS_PACKAGE_NAME "SPAdes 3.16.0-dev")
set(CPACK_NSIS_UNINSTALL_NAME "Uninstall")
set(CPACK_OBJCOPY_EXECUTABLE "/usr/bin/objcopy")
set(CPACK_OBJDUMP_EXECUTABLE "/usr/bin/objdump")
set(CPACK_OUTPUT_CONFIG_FILE "/home/yopican/Downloads/bi/git/splitter/assembler/src/cmake-build-debug/CPackConfig.cmake")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/")
set(CPACK_PACKAGE_DESCRIPTION_FILE "/home/yopican/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/233.13135.93/bin/cmake/linux/x64/share/cmake-3.27/Templates/CPack.GenericDescription.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Project built using CMake")
set(CPACK_PACKAGE_FILE_NAME "SPAdes-3.16.0-dev-Linux")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "SPAdes 3.16.0-dev")
set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "SPAdes 3.16.0-dev")
set(CPACK_PACKAGE_NAME "SPAdes")
set(CPACK_PACKAGE_RELOCATABLE "true")
set(CPACK_PACKAGE_VENDOR "Saint Petersburg State University")
set(CPACK_PACKAGE_VERSION "3.16.0-dev")
set(CPACK_PACKAGE_VERSION_MAJOR "3")
set(CPACK_PACKAGE_VERSION_MINOR "16")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_READELF_EXECUTABLE "/usr/bin/readelf")
set(CPACK_RESOURCE_FILE_LICENSE "/home/yopican/Downloads/bi/git/splitter/assembler/src/../LICENSE")
set(CPACK_RESOURCE_FILE_README "/home/yopican/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/233.13135.93/bin/cmake/linux/x64/share/cmake-3.27/Templates/CPack.GenericDescription.txt")
set(CPACK_RESOURCE_FILE_WELCOME "/home/yopican/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/233.13135.93/bin/cmake/linux/x64/share/cmake-3.27/Templates/CPack.GenericWelcome.txt")
set(CPACK_SET_DESTDIR "OFF")
set(CPACK_SOURCE_GENERATOR "TBZ2")
set(CPACK_SOURCE_IGNORE_FILES "tools;test;web_service;online_vis;cap")
set(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/home/yopican/Downloads/bi/git/splitter/assembler/src/cmake-build-debug/CPackSourceConfig.cmake")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "SPAdes-3.16.0-dev")
set(CPACK_STRIP_FILES "bin/spades-bwa;bin/spades-core;bin/spades-corrector-core;bin/spades-gbuilder;bin/spades-gmapper;bin/spades-hammer;bin/spades-ionhammer;bin/spades-kmercount")
set(CPACK_SYSTEM_NAME "Linux")
set(CPACK_THREADS "1")
set(CPACK_TOPLEVEL_TAG "Linux")
set(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/home/yopican/Downloads/bi/git/splitter/assembler/src/cmake-build-debug/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()
