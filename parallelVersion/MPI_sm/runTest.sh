#!/bin/bash
if [ "$#" -ne 2 ] 
then
  echo "Usage: $0 inputFile compiler"
  exit 1
fi

nloops=3

# Determining MPI implementation and binding options #
MPI=`mpiexec --version | head -1 | awk '{print $1}' ` 

if [ "$MPI" == "HYDRA" ]; then
    echo "MPICH"
    bindings="--bind-to socket"
    export HYDRA_TOPO_DEBUG=1
elif [ "$MPI" == "Intel(R)" ]; then
    echo "Intel MPI"
    bindings="-genv I_MPI_PIN_DOMAIN=core -genv I_MPI_PIN_ORDER=scatter -genv I_MPI_DEBUG=4"
elif [ "$MPI" == "mpiexec" ]; then
    echo "open-mpi"
    bindings="--bind-to core --report-bindings"
fi
# end of Determining MPI implementation and binding options #

npt=`grep -c ^processor /proc/cpuinfo`
numaNodes=`lscpu | grep "NUMA node(s):" | awk '{}{print $3}{}'`
tpc=`lscpu | grep "Thread(s) per core:" | awk '{}{print $4}{}'`
np="$(($npt / $tpc))"
npps="$(($np / $numaNodes))"
npm1="$(($np - 1))"

if [ -n "$PGI" ]; then
    echo "Pgi Compiler"
elif [ -n "$INTEL_LICENSE_FILE" ]; then
    echo "Intel Compiler"
    #np=15
    #npps="$(($np / $numaNodes))"
    #npm1="$(($np - 1))"
else
    echo "Gnu Compiler"
fi


rm -f Mpi_sm_Result.txt
for i in 1 `seq 2 2 $np`; do
    for j in  `seq 1 $nloops`; do
        echo number of processors: $i, run number: $j 
        mpiexec  $bindings -n $i p6 $1 | grep finish >>  Mpi_sm_Result.txt
    done
done

mkdir -p ../../../plots/PureShareMemory/$(hostname)/$2
cat Mpi_sm_Result.txt | awk '{}{print $6, $3}{}' | awk '{Prod[$1]++; min[$1]=Prod[$1]==1||min[$1]>$2?$2:min[$1]} END{ for (var in Prod) printf "%s processors: the min is %f\n", var,min[var]}'  | sort -n   > ../../../plots/PureShareMemory/$(hostname)/$2/p6_MPI_sm.txt

#rm Mpi_sm_Result.txt
#rm RunHistory.dat

