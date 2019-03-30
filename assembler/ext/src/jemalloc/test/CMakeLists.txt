# Most test executables have the same properties, so let's make a common test
# addition function below.
# Also note: some tests need additional libraries for linking or defines set.
# After this function returns, you can inject additional properties to the
# target formatted as: ${PROJECT_NAME}-test-${exename}
function(createTest exename sourcefile testKind)
    set(namespacedTestTarget ${PROJECT_NAME}-test-${exename})

    add_executable(${namespacedTestTarget} ${sourcefile})
    set_target_properties(${namespacedTestTarget} PROPERTIES
        OUTPUT_NAME "${exename}")

    # Depend on jet-static so symbols get generated before tests are built
    add_dependencies(${namespacedTestTarget} ${PROJECT_NAME}-jet-static)

    target_include_directories(${namespacedTestTarget}
        PRIVATE ${PROJECT_SOURCE_DIR}/test/include)

    target_compile_definitions(${namespacedTestTarget}
        PRIVATE ${testKind})

    # Note: we are only testing against the static libararies.
    target_link_libraries(${namespacedTestTarget}
        ${PROJECT_NAME}IntegrationTest
        ${PROJECT_NAME}-static
        ${PROJECT_NAME}-jet-static)

    target_link_libraries(${namespacedTestTarget} pthread dl)

    add_test(NAME ${exename} COMMAND ${namespacedTestTarget})
    add_dependencies(${PROJECT_NAME}-check ${namespacedTestTarget})
endfunction()

add_subdirectory(src)
add_subdirectory(integration)
add_subdirectory(unit)
add_subdirectory(stress)
