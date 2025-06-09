# tests/CMakeLists.txt
# Combine all tests into a single executable and register with CTest

# Enable CTest
enable_testing()

# Single test executable for all tests
add_executable(test_all
  test_GenomicRegionCollection.cpp
  test_BamReader.cpp
)

target_link_libraries(test_all PRIVATE
  seqlib
  Catch2::Catch2
)

# Register the test with CTest
add_test(NAME AllTests COMMAND test_all)

# Ensure coverage instrumentation file is generated when CTest runs
set_tests_properties(AllTests
  PROPERTIES
    ENVIRONMENT "LLVM_PROFILE_FILE=$<TARGET_FILE_DIR:test_all>/tests.profraw"
)
