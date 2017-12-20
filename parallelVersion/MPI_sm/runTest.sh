#!/bin/bash
if [ "$#" -ne 2 ] 
then
  echo "Usage: $0 inputFile compiler"
  exit 1
fi

nloops=5
np=`grep -c ^processor /proc/cpuinfo`


rm -f Mpi_sm_Result.txt
for i in  `seq 1 $np`; do
    for j in  `seq 1 $nloops`; do
        echo number of processors: $i, run number: $j 
        mpiexec -n $i p6 $1 | grep finish >>  Mpi_sm_Result.txt
    done
done

mkdir -p ../../../plots/PureShareMemory/$(hostname)/$2
cat Mpi_sm_Result.txt | awk '{}{print $6, $3}{}' | awk '{Prod[$1]++; min[$1]=Prod[$1]==1||min[$1]>$2?$2:min[$1]} END{ for (var in Prod) printf "%s processors: the min is %f\n", var,min[var]}'  | sort -n   > ../../../plots/PureShareMemory/$(hostname)/$2/p6_MPI_sm.txt

rm Mpi_sm_Result.txt
#rm RunHistory.dat

