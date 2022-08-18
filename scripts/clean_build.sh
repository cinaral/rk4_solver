#!/bin/bash

#* remove build/
rm -r ./build/
if [ $? -eq 0 ]; then
   echo "Removed ./build"
fi

#* remove dat/
rm -r ./dat/ 
if [ $? -eq 0 ]; then
   echo "Removed ./dat"
fi

