cmake_minimum_required(VERSION 3.10)

project(seqlib)

# Set the C++ standard required for the project
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set O2 optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

# Include directories for headers
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/fermi-lite
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa    
)

# Look for htslib on the system
find_package(htslib QUIET)

# Find required system level type libraries
find_package(Threads REQUIRED)  
find_package(ZLIB REQUIRED)

## LZMA
find_path(LZMA_INCLUDE_DIR NAMES lzma.h)
find_library(LZMA_LIBRARY NAMES lzma)

if(NOT LZMA_INCLUDE_DIR OR NOT LZMA_LIBRARY)
  message(FATAL_ERROR "LZMA library or headers not found!")
endif()

## BZip2
find_package(BZip2 REQUIRED)

## HTSLIB
if (htslib_FOUND)
    # If htslib was found on the system, use it
    message(STATUS "Using system htslib")
else()
    set(HTSLIB_DIR "" CACHE PATH "Path to HTSLib root directory")

    if (NOT HTSLIB_DIR)
        message(FATAL_ERROR "HTSLIB_DIR not specified. Please specify -DHTSLIB_DIR=/path/to/htslib")
    else()
        # Automatically set include and library paths based on HTSLIB_DIR
        include_directories(${HTSLIB_DIR}/include)
        link_directories(${HTSLIB_DIR}/lib)
    endif()
endif()    

# Find all source files
file(GLOB SOURCES "src/*.cpp" "src/*.c")

# Silence the deprecation warning in jsoncpp.cpp only
set_source_files_properties(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/jsoncpp.cpp"
  PROPERTIES
    COMPILE_FLAGS "-Wno-deprecated-declarations"
)

# Generate the executable
add_library(seqlib ${SOURCES})

# If the submodules already have Makefiles, you can use custom commands to invoke make
# in those directories.
add_custom_target(
    BuildBWA ALL
    COMMAND make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bwa
)

add_custom_target(
    BuildFermiLite ALL
    COMMAND make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/fermi-lite
)

# If MyProject depends on the submodules, you might want to ensure they're built first
add_dependencies(seqlib BuildBWA BuildFermiLite)

target_link_libraries(seqlib
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/libbwa.a
    ${CMAKE_CURRENT_SOURCE_DIR}/fermi-lite/libfml.a
    Threads::Threads
    ZLIB::ZLIB
    hts
    ${LZMA_LIBRARY}
    BZip2::BZip2
)
