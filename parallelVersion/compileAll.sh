#!/bin/bash


mkdir -p  MPI_sm/buildGnu
mkdir -p  OpenMP/buildGnu

mkdir -p  MPI_sm/buildIntel
mkdir -p  OpenMP/buildIntel

mkdir -p  MPI_sm/buildPgi
mkdir -p  OpenMP/buildPgi



cd OpenMP/buildGnu
cmake .. ; make clean; make
cd ../../MPI_sm/buildGnu
cmake .. ; make clean; make
cd ../../

export CC=icc
export CXX=icpc
source setIcc intel64
source setImpi


cd OpenMP/buildIntel
cmake .. ; make clean; make
cd ../../MPI_sm/buildIntel
cmake .. ; make clean; make
cd ../../

export CC=pgcc
export CXX=pgc++
source setPgi 18.3
source setPgiMpi 18.3

cd OpenMP/buildPgi
cmake .. ; make clean; make
cd ../../MPI_sm/buildPgi
cmake .. ; make clean; make
cd ../../

