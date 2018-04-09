#!/bin/bash
if [ "$#" -ne 2 ]
then
  echo "Usage: $0 inputFile compiler"
  exit 1
fi

nloops=3

npt=`grep -c ^processor /proc/cpuinfo`
numaNodes=`lscpu | grep "NUMA node(s):" | awk '{}{print $3}{}'`
tpc=`lscpu | grep "Thread(s) per core:" | awk '{}{print $4}{}'`
np="$(($npt / $tpc))"
npps="$(($np / $numaNodes))"
npm1="$(($np - 1))"
tps=$((npps*tpc))

seqArray=()
##########################################
k=0
for i in  `seq 0 $((npps-1))`; do
    kk=$((k))
    for j in `seq 0 $((numaNodes-1))`; do
        seqArray[i*$numaNodes+j]=$((kk))
	kk=$((kk+tps))
    done
    k=$((k+tpc))
done
##########################################

#echo ${seqArray[*]}
sequence=''
for p in `seq 0 $((  npm1  ))`; do
    sequence+=${seqArray[p]}','
done
sequence=${sequence%?}
#echo $sequence

if [ -n "$PGI" ]; then
    echo "Pgi Compiler"
elif [ -n "$INTEL_LICENSE_FILE" ]; then
    echo "Intel Compiler"
else
    echo "Gnu Compiler"
fi



rm -f Mpi_sm_Result.txt
for i in 1 `seq 2 2 $np`; do
    for j in  `seq 1 $nloops`; do
        echo number of processors: $i, run number: $j
        aprun -cc $sequence   -n $i -N $i p6 $1 | grep finish >>  Mpi_sm_Result.txt
    done
done

mkdir -p ../../../plots/PureShareMemory/$(hostname)/$2
cat Mpi_sm_Result.txt | awk '{}{print $6, $3}{}' | awk '{Prod[$1]++; min[$1]=Prod[$1]==1||min[$1]>$2?$2:min[$1]} END{ for (var in Prod) printf "%s processors: the min is %f\n", var,min[var]}'  | sort -n   > ../../../plots/PureShareMemory/$(hostname)/$2/p6_MPI_sm.txt

#rm Mpi_sm_Result.txt
#rm RunHistory.dat

