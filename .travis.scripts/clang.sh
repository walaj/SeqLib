#!/bin/bash

## borrowed from Gamgee project: https://github.com/broadinstitute/gamgee/blob/master/.travis_scripts/clang.sh

echo "...building with clang"

wget -O - http://llvm.org/apt/llvm-snapshot.gpg.key | sudo apt-key add -
sudo apt-add-repository 'deb http://llvm.org/apt/precise/ llvm-toolchain-precise-3.5 main'
sudo apt-get -qq update
sudo apt-get -qq --force-yes install clang-3.5 clang-modernize-3.5 # clang-format-3.5
sudo update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-3.5 1
sudo rm /usr/local/clang-3.4/bin/clang++
