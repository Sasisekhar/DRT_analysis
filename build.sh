#!/bin/bash

if [ -d "build" ]; then rm -rf build; fi
mkdir -p build
cd build || exit
rm -rf *
cmake ..
make
cd ..
