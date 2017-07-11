project(read_cloud_statistics CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(read_cloud_statistics main.cpp
                   read_cloud_statistics_extractor.cpp
                    ../spades/repeat_resolving.cpp)

target_link_libraries(read_cloud_statistics ConsensusCore common_modules nlopt BamTools ssw ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
    set_target_properties(read_cloud_statistics PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS read_cloud_statistics
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/read_cloud_statistics"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")