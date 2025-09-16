# Photron - Monte Carlo Photon Transport Renderer
# Modern CMake Build System (CMake 3.20+) - MCML 3.0.0 Strategy
# Author: Enhanced by GitHub Copilot
# Date: 2025
#
# Build Instructions:
#   From build/ # Optional: Create package configuration
set(CPACK_PACKAGE_NAME "Photron")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Monte Carlo Photon Transport Renderer")
set(CPACK_PACKAGE_VENDOR "M.H.J. Lam")
include(CPack)#     cmake -S . -B output
#     cmake --build output --config Release
#   From root directory:
#     cmake -S build -B build/output  
#     cmake --build build/output --config Release
#   For vcpkg manifest mode (if using vcpkg.json):
#     cmake -S . -B output
#     cmake --build output -- /p:VcpkgEnableManifest=true
# This puts all intermediate files in build/output/ and final executable in bin/

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

# Enable vcpkg manifest mode if vcpkg.json is present
if(DEFINED ENV{VCPKG_ROOT})
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")
    # Enable manifest mode and point to the correct vcpkg.json location
    if(EXISTS "${CMAKE_SOURCE_DIR}/vcpkg.json")
        set(VCPKG_MANIFEST_DIR "${CMAKE_SOURCE_DIR}" CACHE STRING "Vcpkg manifest directory")
        set(VCPKG_MANIFEST_MODE ON CACHE BOOL "Enable vcpkg manifest mode")
        # For Visual Studio generator, automatically enable vcpkg manifest in MSBuild
        if(CMAKE_GENERATOR MATCHES "Visual Studio")
            set(CMAKE_VS_GLOBALS "VcpkgEnableManifest=true")
        endif()
    endif()
endif()

# Enable modern CMake policies
cmake_policy(SET CMP0074 NEW)  # find_package uses <PackageName>_ROOT variables

# Project directories - following MCML 3.0.0 strategy
# When using -S . -B output, CMAKE_SOURCE_DIR is build/ and CMAKE_BINARY_DIR is build/output/
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(BIN_DIR ${CMAKE_SOURCE_DIR}/bin)  # Always relative to the source dir (build/)
set(CONFIG_DIR ${CMAKE_SOURCE_DIR}/config)
set(BUILD_OUTPUT_DIR ${CMAKE_BINARY_DIR})  # This will be build/output/

# Create necessary directories
file(MAKE_DIRECTORY ${BIN_DIR})

# Find packages using vcpkg
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};$ENV{VCPKG_ROOT}/installed/x64-windows/share")

find_package(glm CONFIG REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(glew CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)
find_package(cxxopts CONFIG REQUIRED)

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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996 /wd4201")  # Disable deprecation and nameless struct/union warnings
    set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od /DDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG")
endif()

# Create output directories
file(MAKE_DIRECTORY ${BIN_DIR})

# Collect all source files automatically
file(GLOB_RECURSE PHOTRON_SOURCES
    ${SRC_DIR}/main.cpp
    ${SRC_DIR}/app.cpp
    ${SRC_DIR}/math/*.cpp
    ${SRC_DIR}/renderer/*.cpp
    ${SRC_DIR}/simulator/*.cpp
)

# Combine Photron and MCML sources
set(ALL_SOURCES ${PHOTRON_SOURCES})

# Collect all header files (for IDE integration)
file(GLOB_RECURSE PHOTRON_HEADERS
    ${SRC_DIR}/*.hpp
    ${SRC_DIR}/math/*.hpp
    ${SRC_DIR}/renderer/*.hpp
    ${SRC_DIR}/simulator/*.hpp
)

# Print discovered files for debugging
message(STATUS "Discovered Photron source files:")
foreach(src_file ${PHOTRON_SOURCES})
    message(STATUS "  ${src_file}")
endforeach()

# Create executable
add_executable(Photron ${ALL_SOURCES} ${PHOTRON_HEADERS})

# Simple target properties - exe goes to root bin, intermediates to build/output
set_target_properties(Photron PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${BIN_DIR}
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${BIN_DIR}
)

# Add include directories
target_include_directories(Photron PRIVATE 
    ${SRC_DIR}
    ${SRC_DIR}/structs
    ${SRC_DIR}/renderer
    ${SRC_DIR}/simulator
)

# Link libraries using vcpkg targets
target_link_libraries(Photron PRIVATE 
    glfw 
    GLEW::GLEW 
    opengl32 
    glm::glm 
    imgui::imgui
    cxxopts::cxxopts
)

# Platform-specific additional libraries
if(WIN32)
    target_link_libraries(Photron PRIVATE opengl32)
elseif(APPLE)
    target_link_libraries(Photron PRIVATE 
        "-framework OpenGL"
        "-framework Cocoa"
        "-framework IOKit"
        "-framework CoreVideo"
    )
else()
    # Linux and other Unix-like systems
    target_link_libraries(Photron PRIVATE GL GLU)
endif()

# Set output directories following MCML 3.0.0 strategy
set_target_properties(Photron PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${BIN_DIR}
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${BIN_DIR}
    # Intermediate files (object files) go to output directory
    ARCHIVE_OUTPUT_DIRECTORY ${BUILD_OUTPUT_DIR}
    LIBRARY_OUTPUT_DIRECTORY ${BUILD_OUTPUT_DIR}
    # Let vcpkg automatically copy required DLLs
    VS_GLOBAL_VcpkgEnabled true
    VS_GLOBAL_VcpkgEnableManifest true
)

# Enable automatic DLL copying for vcpkg dependencies
if(DEFINED ENV{VCPKG_ROOT})
    set_target_properties(Photron PROPERTIES
        VCPKG_APPLOCAL_DEPS ON
    )
endif()

# Copy configuration files to main bin directory
add_custom_command(TARGET Photron POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CONFIG_DIR}/" "${BIN_DIR}/config/"
)

# Copy shader files to main bin directory
add_custom_command(TARGET Photron POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/shaders/" "${BIN_DIR}/shaders/"
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
install(TARGETS Photron
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
