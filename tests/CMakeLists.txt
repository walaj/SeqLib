# tests/CMakeLists.txt
# create an executable for each test file


add_executable(test_all
  test_GenomicRegionCollection.cpp
  test_BamReader.cpp
)

target_link_libraries(test_all PRIVATE
  seqlib
  Catch2::Catch2
)

# Register the test executable with CTest
add_test(NAME AllTests COMMAND test_all)
