#!/bin/sh
#
# This file needs to be modified according to the computational environment.
#
# WARNING: This test will take a VERY LONG time!
#

MPIRUN=mpirun
NP=-np

rm -f test*.out

make > /dev/null

# 1x1 processor grid
$MPIRUN $NP 1 ./EXRAND1 > test1_1x1.out
$MPIRUN $NP 1 ./EXRAND2 > test2_1x1.out

# 2x2 processor grid
$MPIRUN $NP 4 ./EXRAND1 > test1_2x2.out
$MPIRUN $NP 4 ./EXRAND2 > test2_2x2.out

# 4x4 processor grid
$MPIRUN $NP 16 ./EXRAND1 > test1_4x4.out

grep "ALL OK" *.out >  summary.txt
CNT="`cat summary.txt | wc -l`" 
echo $CNT 'tests of 5 passed'
