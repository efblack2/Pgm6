#!/bin/bash

cd OpenMP/buildGnu
make clean; make
../runTest.sh ../../examples/lowRes.txt gnu
cd ../../MPI_sm/buildGnu
make clean; make
../runTest.sh ../../examples/lowRes.txt gnu

cd ../../

source setIcc intel64
source setImpi

cd OpenMP/buildIntel
make clean; make
../runTest.sh ../../examples/lowRes.txt intel
cd ../../MPI_sm/buildIntel
make clean; make
../runTest.sh ../../examples/lowRes.txt intel

cd ../../

source setPgi 18.3
source setPgiMpi 18.3


cd OpenMP/buildPgi
make clean; make
../runTest.sh ../../examples/lowRes.txt pgi
cd ../../MPI_sm/buildPgi
make clean; make
../runTest.sh ../../examples/lowRes.txt pgi

cd ../../

