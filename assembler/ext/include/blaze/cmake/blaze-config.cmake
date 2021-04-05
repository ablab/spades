#---------------------------------------------------------------------------
#
# blaze-config.cmake - CMake configuration file for external projects.
# Use this by invoking
#
#   find_package(blaze)
#
# The module defines blaze::blaze IMPORTED target

include("${CMAKE_CURRENT_LIST_DIR}/blaze-targets.cmake")
message(STATUS "Found blaze")
