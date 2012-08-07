#!/bin/csh -f

set input = $1
set output = $2
set ref = $3
set procs = $4

setenv IMAGIC_BATCH 1
echo "! "
echo "! "
echo "! ====================== "
echo "! IMAGIC ACCUMULATE FILE "
echo "! ====================== "
echo "! "
echo "! "
echo "! IMAGIC program: mralign ----------------------------------------------"
echo "! "
/opt/qb3/imagic-110326/openmpi/bin/mpirun -np 8 -x IMAGIC_BATCH  /opt/qb3/imagic-110326/align/mralign.e_mpi <<EOF
YES
$procs
FRESH
ALL_REFERENCES
ALIGNMENT
BOTH (ROT AND TRANS)
ROTATION_FIRST
CCF
$input
$output
$input
$ref
NO_FILTER
0.2
-180,180
LOW
0.0,0.7
5
NO
EOF
