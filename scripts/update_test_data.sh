#!/bin/bash

#* specifytest data
test_data_list=(
	'./dat/test-rk4_solver-sine-t_arr.dat' 
	'./dat/test-rk4_solver-sine-x_arr_chk.dat'
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
