cmake ..   -DCMAKE_BUILD_TYPE=Debug   -DCMAKE_C_FLAGS="-g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer"   -DCMAKE_CXX_FLAGS="-g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer"   -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined" -DHTSLIB_DIR=/opt/homebrew/Cellar/htslib/1.21 -DENABLE_COVERAGE=ON

## run tests
LLVM_PROFILE_FILE="tests.profraw" ctest --output-on-failure

## merge raw data
/opt/homebrew/opt/llvm/bin/llvm-profdata merge -sparse tests.profraw -o tests.profdata

# generate HTML
/opt/homebrew/opt/llvm/bin/llvm-cov show ./tests/test_all \
  -instr-profile=tests.profdata \
  -format=html \
  -output-dir=coverage-report \
  -Xdemangler c++filt \
  -Xdemangler -n

