#!/bin/bash

#* project relative path
PROJECT_PATH=..

#*  change the cwd to the script dir temporarily, and hide pushd popd output
pushd () { 
	command pushd "$@" > /dev/null 
}
popd () { 
	command popd "$@" > /dev/null 
}
pushd "$(dirname ${BASH_SOURCE:0})"
trap popd EXIT #*

pushd $PROJECT_PATH/build/bin
./step-benchmark.exe
echo ""
./loop-benchmark.exe
echo ""

echo "$0 done."
