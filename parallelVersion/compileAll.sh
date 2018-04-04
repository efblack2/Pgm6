#!/bin/bash

rm -rf ~/RunHistory.dat
touch ~/RunHistory.dat

mkdir -p  MPI_sm/buildGnu
ln -s ~/RunHistory.dat MPI_sm/buildGnu/RunHistory.dat
mkdir -p  OpenMP/buildGnu
ln -s ~/RunHistory.dat OpenMP/buildGnu/RunHistory.dat

mkdir -p  MPI_sm/buildIntel
ln -s ~/RunHistory.dat MPI_sm/buildIntel/RunHistory.dat
mkdir -p  OpenMP/buildIntel
ln -s ~/RunHistory.dat OpenMP/buildIntel/RunHistory.dat

mkdir -p  MPI_sm/buildPgi
ln -s ~/RunHistory.dat MPI_sm/buildPgi/RunHistory.dat
mkdir -p  OpenMP/buildPgi
ln -s ~/RunHistory.dat OpenMP/buildPgi/RunHistory.dat


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

