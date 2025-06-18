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

# Enhanced SuiteSparse AMD detection
message(STATUS "Searching for SuiteSparse AMD...")

# Method 1: Try pkg-config first
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(AMD QUIET amd)
    if(NOT AMD_FOUND)
        pkg_check_modules(AMD QUIET suitesparse)
    endif()
    if(AMD_FOUND)
        message(STATUS "Found AMD via pkg-config")
    endif()
endif()

# Method 2: Manual search with comprehensive paths
if(NOT AMD_FOUND)
    message(STATUS "pkg-config failed, trying manual search...")
    
    # Common SuiteSparse installation paths
    set(SUITESPARSE_SEARCH_PATHS
        /usr/include/suitesparse
        /usr/local/include/suitesparse
        /opt/local/include/suitesparse
        /usr/include
        /usr/local/include
        /opt/homebrew/include/suitesparse  # macOS Homebrew
        /opt/homebrew/include
        /usr/local/Cellar/suite-sparse/*/include  # macOS Homebrew versioned
        C:/vcpkg/installed/x64-windows/include/suitesparse  # Windows vcpkg
        $ENV{CONDA_PREFIX}/include/suitesparse  # Conda environments
        $ENV{CONDA_PREFIX}/include
    )
    
    set(SUITESPARSE_LIB_PATHS
        /usr/lib
        /usr/local/lib
        /opt/local/lib
        /usr/lib/x86_64-linux-gnu
        /usr/lib64
        /opt/homebrew/lib  # macOS Homebrew
        /usr/local/Cellar/suite-sparse/*/lib  # macOS Homebrew versioned
        C:/vcpkg/installed/x64-windows/lib  # Windows vcpkg
        $ENV{CONDA_PREFIX}/lib  # Conda environments
    )
    
    # Find AMD header
    find_path(AMD_INCLUDE_DIR 
        NAMES amd.h
        PATHS ${SUITESPARSE_SEARCH_PATHS}
        DOC "AMD include directory"
    )
    
    # Find AMD library
    find_library(AMD_LIBRARY 
        NAMES amd libamd
        PATHS ${SUITESPARSE_LIB_PATHS}
        DOC "AMD library"
    )
    
    # Find SuiteSparseConfig library (often required)
    find_library(SUITESPARSECONFIG_LIBRARY 
        NAMES suitesparseconfig libsuitesparseconfig
        PATHS ${SUITESPARSE_LIB_PATHS}
        DOC "SuiteSparseConfig library"
    )
    
    # Also try to find BLAS (sometimes required by SuiteSparse)
    find_library(BLAS_LIBRARY 
        NAMES blas libblas openblas libopenblas
        PATHS ${SUITESPARSE_LIB_PATHS}
        DOC "BLAS library"
    )
    
    # Also try to find LAPACK (sometimes required by SuiteSparse)
    find_library(LAPACK_LIBRARY 
        NAMES lapack liblapack
        PATHS ${SUITESPARSE_LIB_PATHS}
        DOC "LAPACK library"
    )
    
    message(STATUS "AMD search results:")
    message(STATUS "  AMD_INCLUDE_DIR: ${AMD_INCLUDE_DIR}")
    message(STATUS "  AMD_LIBRARY: ${AMD_LIBRARY}")
    message(STATUS "  SUITESPARSECONFIG_LIBRARY: ${SUITESPARSECONFIG_LIBRARY}")
    message(STATUS "  BLAS_LIBRARY: ${BLAS_LIBRARY}")
    message(STATUS "  LAPACK_LIBRARY: ${LAPACK_LIBRARY}")
    
    if(AMD_INCLUDE_DIR AND AMD_LIBRARY)
        set(AMD_FOUND TRUE)
        set(AMD_INCLUDE_DIRS ${AMD_INCLUDE_DIR})
        set(AMD_LIBRARIES ${AMD_LIBRARY})
        
        # Add SuiteSparseConfig if found
        if(SUITESPARSECONFIG_LIBRARY)
            list(APPEND AMD_LIBRARIES ${SUITESPARSECONFIG_LIBRARY})
        endif()
        
        # Add BLAS/LAPACK if found and needed
        if(BLAS_LIBRARY)
            list(APPEND AMD_LIBRARIES ${BLAS_LIBRARY})
        endif()
        if(LAPACK_LIBRARY)
            list(APPEND AMD_LIBRARIES ${LAPACK_LIBRARY})
        endif()
        
        message(STATUS "Found AMD manually:")
        message(STATUS "  Include dir: ${AMD_INCLUDE_DIR}")
        message(STATUS "  Libraries: ${AMD_LIBRARIES}")
    endif()
endif()

# Method 3: Try system package manager hints
if(NOT AMD_FOUND)
    message(STATUS "Manual search failed, checking system package hints...")
    
    # Check if we can find any SuiteSparse files
    execute_process(
        COMMAND find /usr -name "amd.h" 2>/dev/null
        OUTPUT_VARIABLE AMD_FIND_OUTPUT
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    if(AMD_FIND_OUTPUT)
        message(STATUS "Found amd.h at: ${AMD_FIND_OUTPUT}")
        get_filename_component(FOUND_AMD_DIR ${AMD_FIND_OUTPUT} DIRECTORY)
        set(AMD_INCLUDE_DIR ${FOUND_AMD_DIR})
        
        # Try to find corresponding library
        get_filename_component(FOUND_AMD_BASE ${FOUND_AMD_DIR} DIRECTORY)
        find_library(AMD_LIBRARY 
            NAMES amd libamd
            PATHS ${FOUND_AMD_BASE}/lib ${FOUND_AMD_BASE}/lib64
            NO_DEFAULT_PATH
        )
        
        if(AMD_LIBRARY)
            set(AMD_FOUND TRUE)
            set(AMD_INCLUDE_DIRS ${AMD_INCLUDE_DIR})
            set(AMD_LIBRARIES ${AMD_LIBRARY})
            message(STATUS "Found AMD via system search:")
            message(STATUS "  Include dir: ${AMD_INCLUDE_DIR}")
            message(STATUS "  Library: ${AMD_LIBRARY}")
        endif()
    endif()
endif()

# Final status report
if(AMD_FOUND)
    message(STATUS "✓ SuiteSparse AMD found - enabling high-performance minimum degree ordering")
    message(STATUS "  AMD include dirs: ${AMD_INCLUDE_DIRS}")
    message(STATUS "  AMD libraries: ${AMD_LIBRARIES}")
else()
    message(STATUS "✗ SuiteSparse AMD not found - will use fallback implementation")
    message(STATUS "")
    message(STATUS "To install SuiteSparse AMD:")
    message(STATUS "  Ubuntu/Debian: sudo apt-get install libsuitesparse-dev")
    message(STATUS "  CentOS/RHEL:   sudo yum install suitesparse-devel")
    message(STATUS "  macOS:         brew install suite-sparse")
    message(STATUS "  conda:         conda install -c conda-forge suitesparse")
    message(STATUS "")
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
    target_compile_definitions(hypergraph_ordering_lib PUBLIC HAVE_AMD)
    target_include_directories(hypergraph_ordering_lib PUBLIC ${AMD_INCLUDE_DIRS})
    target_link_libraries(hypergraph_ordering_lib PUBLIC ${AMD_LIBRARIES})
    
    # Test if we can compile with AMD
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_INCLUDES ${AMD_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${AMD_LIBRARIES})
    
    # Save current flags
    set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -I${AMD_INCLUDE_DIRS}")
    
    check_cxx_source_compiles("
        extern \"C\" {
        #include <amd.h>
        }
        int main() {
            double Control[AMD_CONTROL], Info[AMD_INFO];
            amd_defaults(Control);
            return 0;
        }
    " AMD_COMPILES)
    
    # Restore flags
    set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
    
    if(NOT AMD_COMPILES)
        message(WARNING "AMD found but cannot compile test program. May need additional libraries.")
        message(STATUS "Consider installing additional dependencies:")
        message(STATUS "  Ubuntu/Debian: sudo apt-get install libblas-dev liblapack-dev")
        message(STATUS "  CentOS/RHEL:   sudo yum install blas-devel lapack-devel")
        
        # Try alternative include paths
        message(STATUS "Trying alternative include paths...")
        find_path(AMD_INCLUDE_DIR_ALT
            NAMES amd.h
            PATHS 
                /usr/include/suitesparse
                /usr/local/include/suitesparse
                /opt/local/include/suitesparse
                /usr/include
                /usr/local/include
            NO_DEFAULT_PATH
        )
        
        if(AMD_INCLUDE_DIR_ALT AND NOT "${AMD_INCLUDE_DIR_ALT}" STREQUAL "${AMD_INCLUDE_DIR}")
            message(STATUS "Found alternative AMD include path: ${AMD_INCLUDE_DIR_ALT}")
            set(AMD_INCLUDE_DIRS ${AMD_INCLUDE_DIR_ALT})
            target_include_directories(hypergraph_ordering_lib PUBLIC ${AMD_INCLUDE_DIRS})
        endif()
    else()
        message(STATUS "✓ AMD compilation test passed")
    endif()
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
    COMMAND ${CMAKE_COMMAND} -E echo "AMD found: ${AMD_FOUND}"
    COMMAND ${CMAKE_COMMAND} -E echo "AMD include dirs: ${AMD_INCLUDE_DIRS}"
    COMMAND ${CMAKE_COMMAND} -E echo "AMD libraries: ${AMD_LIBRARIES}"
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
message(STATUS "")
message(STATUS "=== Build Configuration Summary ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "KaHyPar version: ${KAHYPAR_VERSION}")
message(STATUS "OpenMP support: ${OpenMP_FOUND}")
message(STATUS "AMD support: ${AMD_FOUND}")
message(STATUS "")

# Example usage commands (printed during configure)
message(STATUS "Build commands:")
message(STATUS "  mkdir build && cd build")
message(STATUS "  cmake ..")
message(STATUS "  make -j$(nproc)")
message(STATUS "")
message(STATUS "Run:")
message(STATUS "  ./hypergraph_ordering <matrix.mtx>")
message(STATUS "")