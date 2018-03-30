#!/bin/bash
if [ "$#" -lt 1 ]
then
  echo "Usage: $0  computer"
  exit 1
fi

cd $1/gnu
paste p6_MPI_sm.txt  p6_OpenMP.txt > results.txt

cd ../intel
paste p6_MPI_sm.txt  p6_OpenMP.txt > results.txt

cd ../pgi
paste p6_MPI_sm.txt  p6_OpenMP.txt > results.txt

cd ../..


gnuplot -c plot.gnp $1 
gnuplot -c plotRatio.gnp $1 

mv $1.pdf temp.pdf
pdfunite temp.pdf $1Ratio.pdf $1.pdf

rm temp.pdf $1Ratio.pdf

gnuplot -c plotGnu.gnp $1  
gnuplot -c plotGnuRatio.gnp $1  


mv $1_Gnu.pdf temp.pdf
pdfunite temp.pdf $1Ratio_Gnu.pdf $1_Gnu.pdf

rm temp.pdf $1Ratio_Gnu.pdf




rm `find . -name  results.txt`
