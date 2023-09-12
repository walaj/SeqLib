cmake_minimum_required(VERSION 3.10)

project(seqlib)

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
    # Option for custom include directory
    set(HTSLIB_INCLUDE_DIR "" CACHE PATH "Path to custom htslib include directory")

    # Option for custom library
    set(HTSLIB_LIBRARY "" CACHE FILEPATH "Path to custom htslib library file")

    # Check if the user provided custom paths and they are valid
    if (HTSLIB_INCLUDE_DIR AND EXISTS "${HTSLIB_INCLUDE_DIR}" AND HTSLIB_LIBRARY AND EXISTS "${HTSLIB_LIBRARY}")
        include_directories(${HTSLIB_INCLUDE_DIR})
        message(STATUS "Using custom htslib from specified paths")
    else()
        message(STATUS "htslib not found on system. You can specify custom paths with -DHTSLIB_INCLUDE_DIR and -DHTSLIB_LIBRARY")
    endif()
endif()


# Find all source files
file(GLOB SOURCES "src/*.cpp" "src/*.c")

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
    ${HTSLIB_LIBRARY}
    ${LZMA_LIBRARY}
    BZip2::BZip2
)
