#!/bin/bash

#* project relative path
PROJECT_PATH=../..
CLEAR_DIR=$PROJECT_PATH/build

#*  change the cwd to the script dir temporarily, and hide pushd popd output
pushd () { 
	command pushd "$@" > /dev/null 
}
popd () { 
	command popd "$@" > /dev/null 
}
pushd "$(dirname ${BASH_SOURCE:0})"
trap popd EXIT #*

#* remove dir
if [[ -d "$CLEAR_DIR" ]]; then
rm -r -f $CLEAR_DIR

if [[ $? -eq 0 ]]; then
   echo "Removed $CLEAR_DIR"
fi
fi

echo "$0 done."
