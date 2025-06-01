cmake_minimum_required(VERSION 2.8.12)

# Set policy for project VERSION
if(POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif()

# EIS Framework as part of lmfit project
project(eis_framework LANGUAGES C)

# Inherit C standard from parent
if(NOT CMAKE_C_STANDARD)
    set(CMAKE_C_STANDARD 99)
endif()

# Options
option(BUILD_EIS_EXAMPLES "Build EIS example programs" ON)

# Include directories
include_directories(src)
include_directories(${CMAKE_SOURCE_DIR}/lib)

# Core library sources
set(EIS_CORE_SOURCES
    src/eis_core.c
    src/eis_fitting.c
    src/eis_output.c
)

# Model sources
set(EIS_MODEL_SOURCES
    models/randles_model.c
    models/battery_model.c
)

# Create EIS library
add_library(eis_framework ${EIS_CORE_SOURCES} ${EIS_MODEL_SOURCES})

# Link with lmfit library from parent project (corrected name)
target_link_libraries(eis_framework lmfit)

# Math library (not needed on Windows with MSVC)
if(NOT WIN32 OR NOT MSVC)
    target_link_libraries(eis_framework m)
endif()

# Examples
if(BUILD_EIS_EXAMPLES)
    # Main example
    add_executable(eis_main_example examples/main_example.c)
    target_link_libraries(eis_main_example eis_framework)
    
    # Randles demo
    add_executable(eis_randles_demo examples/randles_demo.c)
    target_link_libraries(eis_randles_demo eis_framework)
    
    # Battery demo (original eis.c model)
    add_executable(eis_battery_demo examples/battery_demo.c)
    target_link_libraries(eis_battery_demo eis_framework)
    
    # Set output directory to be the same as other executables
    if(MSVC)
        set_target_properties(eis_main_example eis_randles_demo eis_battery_demo
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
        )
    endif()
endif()

# Installation (only if LIB_INSTALL is ON)
if(LIB_INSTALL)
    install(TARGETS eis_framework
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
    )
    
    install(FILES src/eis_framework.h
        DESTINATION include
    )
endif()

# Print configuration
message(STATUS "EIS Framework Configuration:")
message(STATUS "  Build examples: ${BUILD_EIS_EXAMPLES}")
message(STATUS "  Will link with: lmfit (from parent project)")