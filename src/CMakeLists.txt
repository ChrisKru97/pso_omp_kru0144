add_executable(main main.c)
add_library(pso pso.c)

set_target_properties(
    main

    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build
)

set_target_properties(
    pso

    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build
)

target_link_libraries(
    main PUBLIC
        ${PROJECT_SOURCE_DIR}/src/pso.c
)

target_include_directories(
    main PUBLIC
        ${PROJECT_SOURCE_DIR}/lib/
)
