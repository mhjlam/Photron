# Photron - Monte Carlo Photon Transport Renderer
# Modern CMake Build System (CMake 3.20+)
# Author: Enhanced by GitHub Copilot
# Date: 2025

cmake_minimum_required(VERSION 3.20)
project(Photron VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Enable CMake integration with vcpkg
if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")
endif()

# Enable modern CMake policies
cmake_policy(SET CMP0074 NEW)  # find_package uses <PackageName>_ROOT variables

# Project directories (relative to build directory)
set(SRC_DIR ${CMAKE_SOURCE_DIR}/../src)
set(BIN_DIR ${CMAKE_SOURCE_DIR}/../bin)
set(CONFIG_DIR ${CMAKE_SOURCE_DIR}/../config)
set(OUTPUT_DIR ${CMAKE_BINARY_DIR}/output)

# Configure build system to use output directory for ALL intermediate files
set(CMAKE_BINARY_DIR ${OUTPUT_DIR})

# Configure build output directories to keep build folder clean
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTPUT_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${OUTPUT_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}/bin)

# For multi-config generators (Visual Studio) - put ALL files in output
foreach(CONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG} CONFIG_UPPER)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG_UPPER} ${OUTPUT_DIR}/${CONFIG}/lib)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG_UPPER} ${OUTPUT_DIR}/${CONFIG}/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG_UPPER} ${OUTPUT_DIR}/${CONFIG}/bin)
endforeach()

# Force Visual Studio to place project files and intermediate files in output directory
if(CMAKE_GENERATOR MATCHES "Visual Studio")
    set(CMAKE_VS_INCLUDE_INSTALL_TO_DEFAULT_BUILD FALSE)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

# Create output directories
file(MAKE_DIRECTORY ${OUTPUT_DIR})
file(MAKE_DIRECTORY ${OUTPUT_DIR}/bin)
file(MAKE_DIRECTORY ${OUTPUT_DIR}/lib)
file(MAKE_DIRECTORY ${BIN_DIR})

# Find packages using vcpkg
find_package(OpenGL REQUIRED)
find_package(glfw3 CONFIG)
find_package(GLEW CONFIG)
find_package(glm CONFIG)

# Compiler-specific flags
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wconversion -Wsign-conversion")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wconversion -Wsign-conversion")
    # Suppress specific warnings that are too noisy
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-utf8 -Wno-deprecated-declarations")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /EHsc")
    # Disable specific MSVC warnings that are too noisy
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996")  # Disable deprecation warnings
    set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od /DDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG")
endif()

# Create output directories
file(MAKE_DIRECTORY ${BIN_DIR})

# Collect all source files automatically
file(GLOB_RECURSE PHOTRON_SOURCES
    ${SRC_DIR}/main.cpp
    ${SRC_DIR}/renderer/*.cpp
    ${SRC_DIR}/simulator/initializer.cpp
    ${SRC_DIR}/simulator/mcml_integration.cpp
    ${SRC_DIR}/simulator/simulator.cpp
    ${SRC_DIR}/utilities/*.cpp
)

# Collect MCML sources (to be compiled in separate namespace)
file(GLOB_RECURSE MCML_SOURCES
    ${SRC_DIR}/mcml/*.cpp
)

# Remove old/backup files and exclude MCML directory from main sources
list(FILTER PHOTRON_SOURCES EXCLUDE REGEX ".*_old\\.cpp$")
list(FILTER PHOTRON_SOURCES EXCLUDE REGEX ".*/mcml/.*")
# Use only the main renderer, exclude all alternative versions
list(FILTER PHOTRON_SOURCES EXCLUDE REGEX ".*/renderer_modern\\.cpp$")
list(FILTER PHOTRON_SOURCES EXCLUDE REGEX ".*/renderer_working\\.cpp$")

# Remove main.cpp from MCML sources to avoid duplicate main functions
list(FILTER MCML_SOURCES EXCLUDE REGEX ".*/main\\.cpp$")

# Combine Photron and MCML sources
set(ALL_SOURCES ${PHOTRON_SOURCES} ${MCML_SOURCES})

# Collect all header files (for IDE integration)
file(GLOB_RECURSE PHOTRON_HEADERS
    ${SRC_DIR}/structs/*.hpp
    ${SRC_DIR}/renderer/*.hpp
    ${SRC_DIR}/simulator/*.hpp
    ${SRC_DIR}/utilities/*.hpp
    ${SRC_DIR}/render/*.hpp
    ${SRC_DIR}/mcml/*.hpp
)

# Print discovered files for debugging
message(STATUS "Discovered Photron source files:")
foreach(src_file ${PHOTRON_SOURCES})
    message(STATUS "  ${src_file}")
endforeach()

message(STATUS "Discovered MCML source files:")
foreach(src_file ${MCML_SOURCES})
    message(STATUS "  ${src_file}")
endforeach()

# Create executable
add_executable(photron ${ALL_SOURCES} ${PHOTRON_HEADERS})

# Set target-specific output directories to organize build files
set_target_properties(photron PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}/bin
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${OUTPUT_DIR}/Release/bin
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${OUTPUT_DIR}/Debug/bin
    ARCHIVE_OUTPUT_DIRECTORY ${OUTPUT_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${OUTPUT_DIR}/Release/lib
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${OUTPUT_DIR}/Debug/lib
)

# Add include directories
target_include_directories(photron PRIVATE 
    ${SRC_DIR}
    ${SRC_DIR}/structs
    ${SRC_DIR}/renderer
    ${SRC_DIR}/simulator
    ${SRC_DIR}/utilities
    ${SRC_DIR}/render
    ${SRC_DIR}/mcml
)

# Link libraries using vcpkg targets
set(PHOTRON_LIBRARIES ${OPENGL_LIBRARIES})

# Try to link modern OpenGL libraries if available
if(TARGET glfw)
    list(APPEND PHOTRON_LIBRARIES glfw)
    message(STATUS "Using GLFW from vcpkg")
endif()

if(TARGET GLEW::GLEW)
    list(APPEND PHOTRON_LIBRARIES GLEW::GLEW)
    message(STATUS "Using GLEW from vcpkg")
endif()

if(TARGET glm::glm)
    list(APPEND PHOTRON_LIBRARIES glm::glm)
    message(STATUS "Using GLM from vcpkg")
endif()

target_link_libraries(photron PRIVATE ${PHOTRON_LIBRARIES})

# Platform-specific additional libraries
if(WIN32)
    target_link_libraries(photron PRIVATE 
        opengl32
    )
elseif(APPLE)
    target_link_libraries(photron PRIVATE 
        "-framework OpenGL"
        "-framework Cocoa"
        "-framework IOKit"
        "-framework CoreVideo"
    )
else()
    # Linux and other Unix-like systems
    target_link_libraries(photron PRIVATE 
        GL 
        GLU
    )
endif()

# Set output directories and intermediate file locations
set_target_properties(photron PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG   "${BIN_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${BIN_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${BIN_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL "${BIN_DIR}"
    # Force Visual Studio to put intermediate files in output directory
    COMPILE_PDB_OUTPUT_DIRECTORY_DEBUG "${OUTPUT_DIR}/Debug"
    COMPILE_PDB_OUTPUT_DIRECTORY_RELEASE "${OUTPUT_DIR}/Release"
    VS_INTERMEDIATE_DIRECTORY_DEBUG "${OUTPUT_DIR}/Debug"
    VS_INTERMEDIATE_DIRECTORY_RELEASE "${OUTPUT_DIR}/Release"
    # Object file directory
    OBJECT_DEPENDS ""
)

# Copy final executable from build/output to main bin directory for easy access
add_custom_command(TARGET photron POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:photron>" "${BIN_DIR}/photron.exe"
    COMMENT "Copying executable from build/output to main bin directory"
)

# Copy configuration files to main bin directory
add_custom_command(TARGET photron POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CONFIG_DIR}/" "${BIN_DIR}/config/"
    COMMENT "Copying configuration files to target directory"
)

# Copy shader files to main bin directory
add_custom_command(TARGET photron POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/../shaders/" "${BIN_DIR}/shaders/"
    COMMENT "Copying shader files to target directory"
)

# Print build information
message(STATUS "Photron v${PROJECT_VERSION} Build Configuration:")
message(STATUS "  Platform:     ${CMAKE_SYSTEM_NAME}")
message(STATUS "  Compiler:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Build Type:   ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "  Source Dir:   ${SRC_DIR}")
message(STATUS "  Binary Dir:   ${BIN_DIR}")
message(STATUS "  Build Dir:    ${CMAKE_BINARY_DIR}")
message(STATUS "  Flags:        ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")

# Optional: Create install targets
install(TARGETS photron
    RUNTIME DESTINATION bin
    COMPONENT runtime
)

install(DIRECTORY ${CONFIG_DIR}/
    DESTINATION bin/config
    COMPONENT configuration
    FILES_MATCHING PATTERN "*.in"
)

# Optional: Create package configuration
set(CPACK_PACKAGE_NAME "Photron")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Monte Carlo Photon Transport Renderer")
set(CPACK_PACKAGE_VENDOR "M.H.J. Lam")
include(CPack)
