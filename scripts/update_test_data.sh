#!/bin/bash

#* specifytest data
test_data_list=(
	"test-rk4_solver-motor-u_arr.dat"
	"test-rk4_solver-motor-x_arr_chk.dat"
	"test-rk4_solver-ball-x_arr_chk.dat"
	)

#* remove existing test data
rm -r ../test/dat 
if [ $? -eq 0 ]; then
   echo "Removed ./test/dat"
fi

#* move test data from ./dat to ./test/dat
mkdir ../test/dat

cp ${test_data_list[*]} ../test/dat
if [ $? -eq 0 ]; then
   echo "Updated ./test/dat"
fi
