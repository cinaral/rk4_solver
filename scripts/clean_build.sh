#!/bin/bash

#* remove build/
mv ./build/_deps ./_deps/
rm -r ./build/
mkdir build
mv ./_deps ./build/_deps
if [ $? -eq 0 ]; then
   echo "Removed build/ except build/_deps"
fi

#* remove dat/
rm -r ./dat/ 
if [ $? -eq 0 ]; then
   echo "Removed ./dat"
fi

