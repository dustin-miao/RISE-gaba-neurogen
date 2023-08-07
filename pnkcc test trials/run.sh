#!/bin/sh

cd "$1"
echo "Running simulation in $1"
make clean
rm mySim
rm -R results
mkdir results
make
pwd
./mySim
echo "Running simulation in $1"
cd "/Users/dustinmiao/Documents/school/rise/RISE-gaba-neurogen"