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


seqArray=()
##########################################
for i in  `seq 0 $((npps-1))`; do
    for j in `seq 0 $((numaNodes-1))`; do
        seqArray[i*$numaNodes+j]=$((i+j*npps))
    done
done
##########################################
#for i in `seq 0 $((npm1))`; do
#    seqArray[i]=$i
#done
##########################################
#for i in `seq 0 2 $((npm1))`; do
#    sequence+=$i','
#done
#for i in `seq 1 2 $((npm1))`; do
#    sequence+=$i','
#done
##########################################

#echo ${seqArray[*]}
sequence=''
for p in `seq 0 $((  npm1  ))`; do
    sequence+=${seqArray[p]}','
done
sequence=${sequence%?}
export OMP_DISPLAY_ENV=true
if [ -n "$PGI" ]; then
    echo "Pgi Compiler"
    export MP_BLIST=$sequence
    export MP_BIND="yes"
    #export MP_BLIST="0-$npm1"
    echo $MP_BLIST
elif [ -n "$INTEL_LICENSE_FILE" ]; then
    echo "Intel Compiler"
    #np=15
    #npps="$(($np / $numaNodes))"
    #npm1="$(($np - 1))"
    export OMP_PROC_BIND=spread
    export OMP_PLACES=cores
    #export KMP_AFFINITY=scatter
    # needed to use dissabled in Blue waters
    #export KMP_AFFINITY=disabled
else
    echo "Gnu Compiler"
    export OMP_PROC_BIND=spread
    export OMP_PLACES=sockets
fi

rm -f openMpResult.txt
for i in 1 `seq 2 2 $np`; do
    export OMP_NUM_THREADS=$i
    for j in  `seq 1 $nloops`; do
        echo number of threads: $i, run number: $j
        p6 $1 | grep finish >>  openMpResult.txt
    done
done

mkdir -p ../../../plots/PureShareMemory/$(hostname)/$2
cat openMpResult.txt | awk '{}{print $6, $3}{}' | awk '{Prod[$1]++; min[$1]=Prod[$1]==1||min[$1]>$2?$2:min[$1]} END{ for (var in Prod) printf "%s threads: the min is %f\n", var,min[var]}'  | sort -n  > ../../../plots/PureShareMemory/$(hostname)/$2/p6_OpenMP.txt

rm openMpResult.txt
#rm RunHistory.dat

