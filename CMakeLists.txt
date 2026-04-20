cmake_minimum_required(VERSION 3.14)
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
#
# Detection order (same as the parent svaba build):
#   1. find_package(htslib)    — only works when htslib ships a
#      htslibConfig.cmake (conda/vcpkg).
#   2. pkg-config --exists htslib
#      — system installs from brew / apt ship htslib.pc, which works.
#   3. find_path/find_library on default system search paths.
#   4. Explicit -DHTSLIB_DIR=/path for non-standard installs.
#
# This duplicates the logic in the parent svaba CMakeLists intentionally
# because SeqLib is reused standalone by other projects (plus parent
# add_subdirectory(SeqLib) runs before the parent's own htslib detection
# kicks in, so we can't rely on the parent's values here).
find_package(htslib QUIET)

set(SEQLIB_HTSLIB_FOUND FALSE)
set(SEQLIB_HTSLIB_LIBRARIES "")

if (htslib_FOUND)
    message(STATUS "SeqLib: htslib found via find_package()")
    set(SEQLIB_HTSLIB_FOUND TRUE)
    set(SEQLIB_HTSLIB_LIBRARIES hts)
endif()

if (NOT SEQLIB_HTSLIB_FOUND)
    find_package(PkgConfig QUIET)
    if (PkgConfig_FOUND)
        pkg_check_modules(HTSLIB_PC QUIET htslib)
        if (HTSLIB_PC_FOUND)
            message(STATUS "SeqLib: htslib found via pkg-config (version ${HTSLIB_PC_VERSION})")
            include_directories(${HTSLIB_PC_INCLUDE_DIRS})
            link_directories(${HTSLIB_PC_LIBRARY_DIRS})
            set(SEQLIB_HTSLIB_FOUND TRUE)
            set(SEQLIB_HTSLIB_LIBRARIES ${HTSLIB_PC_LIBRARIES})
        endif()
    endif()
endif()

if (NOT SEQLIB_HTSLIB_FOUND)
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/hts.h
              DOC "Directory containing htslib/hts.h")
    find_library(HTSLIB_LIBRARY NAMES hts
                 DOC "Path to the hts library")
    if (HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
        message(STATUS "SeqLib: htslib headers at ${HTSLIB_INCLUDE_DIR}, lib at ${HTSLIB_LIBRARY}")
        include_directories(${HTSLIB_INCLUDE_DIR})
        set(SEQLIB_HTSLIB_FOUND TRUE)
        set(SEQLIB_HTSLIB_LIBRARIES ${HTSLIB_LIBRARY})
    endif()
endif()

if (NOT SEQLIB_HTSLIB_FOUND)
    set(HTSLIB_DIR "" CACHE PATH "Path to HTSLib root directory")
    if (HTSLIB_DIR)
        message(STATUS "SeqLib: using explicit HTSLIB_DIR=${HTSLIB_DIR}")
        include_directories(${HTSLIB_DIR}/include)
        link_directories(${HTSLIB_DIR}/lib)
        set(SEQLIB_HTSLIB_FOUND TRUE)
        set(SEQLIB_HTSLIB_LIBRARIES hts)
    else()
        message(FATAL_ERROR
            "SeqLib: htslib not found. Install htslib system-wide (e.g. "
            "`brew install htslib` or `apt install libhts-dev`) so "
            "pkg-config / find_library locate it automatically, or pass "
            "-DHTSLIB_DIR=/path/to/htslib to point at a manual build.")
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
    ${SEQLIB_HTSLIB_LIBRARIES}
    ${LZMA_LIBRARY}
    BZip2::BZip2
)

# enable CTest
enable_testing()

# pull in Catch2
include(FetchContent)
FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.x   # or v3.x
)
FetchContent_MakeAvailable(catch2)

add_definitions(-DTEST_DATA_DIR=\"${CMAKE_SOURCE_DIR}/tests/data\")

# add your tests directory
add_subdirectory(tests)



if(ENABLE_COVERAGE)
  message(STATUS "Building with Clang coverage")
  add_compile_options(-fprofile-instr-generate -fcoverage-mapping -g -O0)
  add_link_options(-fprofile-instr-generate -fcoverage-mapping)
endif()
