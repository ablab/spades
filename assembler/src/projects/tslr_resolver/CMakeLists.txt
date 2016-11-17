project(tslr_resolver CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(tslr_resolver main.cpp
        ../spades/pair_info_count.cpp
        ../spades/mismatch_correction.cpp
        ../spades/second_phase_setup.cpp
        ../spades/gap_closer.cpp
        ../spades/distance_estimation.cpp
        ../spades/hybrid_aligning.cpp)

target_include_directories(tslr_resolver PRIVATE ${EXT_DIR}/include/ConsensusCore)
target_link_libraries(tslr_resolver ConsensusCore common_modules nlopt BamTools ssw ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
    set_target_properties(tslr_resolver PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS tslr_resolver
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/tslr_resolver"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")