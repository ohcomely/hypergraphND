cmake_minimum_required(VERSION 3.16)
project(hypergraph_ordering VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler-specific optimizations
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -pedantic")
endif()

# Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Find required packages
find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)

# Find KaHyPar using pkg-config
pkg_check_modules(KAHYPAR REQUIRED IMPORTED_TARGET libkahypar)

# Find SuiteSparse AMD
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(AMD QUIET suitesparseconfig amd)
endif()

# Alternative: Try to find AMD manually if pkg-config fails
if(NOT AMD_FOUND)
    find_path(AMD_INCLUDE_DIR 
        NAMES amd.h
        PATHS 
            /usr/include/suitesparse 
            /usr/local/include/suitesparse
            /opt/local/include/suitesparse 
            /usr/include 
            /usr/local/include
        PATH_SUFFIXES suitesparse
    )
    
    find_library(AMD_LIBRARY 
        NAMES amd
        PATHS 
            /usr/lib 
            /usr/local/lib 
            /opt/local/lib
            /usr/lib/x86_64-linux-gnu
        PATH_SUFFIXES suitesparse
    )
    
    find_library(SUITESPARSECONFIG_LIBRARY 
        NAMES suitesparseconfig
        PATHS 
            /usr/lib 
            /usr/local/lib 
            /opt/local/lib
            /usr/lib/x86_64-linux-gnu
        PATH_SUFFIXES suitesparse
    )
    
    if(AMD_INCLUDE_DIR AND AMD_LIBRARY)
        set(AMD_FOUND TRUE)
        set(AMD_INCLUDE_DIRS ${AMD_INCLUDE_DIR})
        set(AMD_LIBRARIES ${AMD_LIBRARY})
        if(SUITESPARSECONFIG_LIBRARY)
            list(APPEND AMD_LIBRARIES ${SUITESPARSECONFIG_LIBRARY})
        endif()
        message(STATUS "Found AMD manually:")
        message(STATUS "  Include dir: ${AMD_INCLUDE_DIR}")
        message(STATUS "  Libraries: ${AMD_LIBRARIES}")
    else()
        message(STATUS "AMD search results:")
        message(STATUS "  Include dir: ${AMD_INCLUDE_DIR}")
        message(STATUS "  AMD library: ${AMD_LIBRARY}")
        message(STATUS "  SuiteSparseConfig library: ${SUITESPARSECONFIG_LIBRARY}")
    endif()
endif()

if(AMD_FOUND)
    message(STATUS "SuiteSparse AMD found - enabling high-performance minimum degree ordering")
else()
    message(STATUS "SuiteSparse AMD not found - will use fallback implementation")
endif()

# Include directories
include_directories(include)

# Source files
set(SOURCES
    src/sparse_matrix.cpp
    src/hypergraph_ordering.cpp
    src/hypergraph_construction.cpp
    src/hypergraph_decode.cpp
    src/partition.cpp
    src/minimum_degree.cpp
)

# Create library
add_library(hypergraph_ordering_lib STATIC ${SOURCES})

# Link libraries to the library
target_link_libraries(hypergraph_ordering_lib 
    PUBLIC 
        OpenMP::OpenMP_CXX
        PkgConfig::KAHYPAR
)

# Add AMD if available
if(AMD_FOUND)
    target_compile_definitions(hypergraph_ordering_lib PRIVATE HAVE_AMD)
    target_include_directories(hypergraph_ordering_lib PRIVATE ${AMD_INCLUDE_DIRS})
    target_link_libraries(hypergraph_ordering_lib PRIVATE ${AMD_LIBRARIES})
endif()

# Set include directories for the library
target_include_directories(hypergraph_ordering_lib 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Create executable
add_executable(hypergraph_ordering src/main.cpp)

# Link the library to the executable
target_link_libraries(hypergraph_ordering 
    PRIVATE 
        hypergraph_ordering_lib
)

# Compiler-specific warnings and optimizations
target_compile_options(hypergraph_ordering_lib PRIVATE
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic -Wno-unused-parameter>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic -Wno-unused-parameter>
)

# Add custom target for debugging info
add_custom_target(debug_info
    COMMAND ${CMAKE_COMMAND} -E echo "Build type: ${CMAKE_BUILD_TYPE}"
    COMMAND ${CMAKE_COMMAND} -E echo "C++ flags: ${CMAKE_CXX_FLAGS}"
    COMMAND ${CMAKE_COMMAND} -E echo "C++ flags release: ${CMAKE_CXX_FLAGS_RELEASE}"
    COMMAND ${CMAKE_COMMAND} -E echo "KaHyPar found: ${KAHYPAR_FOUND}"
    COMMAND ${CMAKE_COMMAND} -E echo "KaHyPar version: ${KAHYPAR_VERSION}"
    COMMAND ${CMAKE_COMMAND} -E echo "OpenMP found: ${OpenMP_FOUND}"
)

# Installation (optional)
install(TARGETS hypergraph_ordering
    RUNTIME DESTINATION bin
)

install(TARGETS hypergraph_ordering_lib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)

# Print configuration info
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "KaHyPar version: ${KAHYPAR_VERSION}")
message(STATUS "OpenMP support: ${OpenMP_FOUND}")

# Example usage commands (printed during configure)
message(STATUS "")
message(STATUS "Build commands:")
message(STATUS "  mkdir build && cd build")
message(STATUS "  cmake ..")
message(STATUS "  make -j$(nproc)")
message(STATUS "")
message(STATUS "Run:")
message(STATUS "  ./hypergraph_ordering <matrix.mtx>")
message(STATUS "")