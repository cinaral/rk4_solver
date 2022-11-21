#!/bin/bash

#* project relative path
PROJECT_PATH=../..

#* C/CXX compiler absolute path for MSYS2 on Windows
C_COMPILER_PATH="/c/msys64/mingw64/bin/x86_64-w64-mingw32-gcc.exe"
CXX_COMPILER_PATH="/c/msys64/mingw64/bin/x86_64-w64-mingw32-g++.exe"

#*  change the cwd to the script dir temporarily, and hide pushd popd output
pushd () { 
	command pushd "$@" > /dev/null 
}
popd () { 
	command popd "$@" > /dev/null 
}
pushd "$(dirname ${BASH_SOURCE:0})"
trap popd EXIT #*

UNAME=$(uname)

if [[ $UNAME == "Linux" ]] ; then
	cmake -S $PROJECT_PATH/ -B $PROJECT_PATH/build
elif [[ $UNAME == "MSYS"* ]] ; then
	cmake -S $PROJECT_PATH/ -B $PROJECT_PATH/build -DCMAKE_EXPORT_COMPILE_COMMANDS=1 --no-warn-unused-cli -DCMAKE_C_COMPILER:FILEPATH=$C_COMPILER_PATH -DCMAKE_CXX_COMPILER:FILEPATH=$CXX_COMPILER_PATH -G "MinGW Makefiles"
fi

echo "$0 done."
