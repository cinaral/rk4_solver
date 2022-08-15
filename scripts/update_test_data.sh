#!/bin/bash

#* specifytest data
test_data_list=(
	)

#* remove existing test data
rm -r ./test/dat 
if [ $? -eq 0 ]; then
   echo "Removed ./test/dat"
fi

#* move test data from ./dat to ./test/dat
mkdir ./test/dat

cp ${test_data_list[*]} ./test/dat
if [ $? -eq 0 ]; then
   echo "Updated ./test/dat"
fi
