############################################################################
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(logan-extract CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(logan-extract main.cpp)

target_link_libraries(logan-extract common_modules toolchain ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(logan-extract PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS logan-extract
        DESTINATION bin
        COMPONENT tools)
