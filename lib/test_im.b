#!/bin/csh -f

set tot = $1
set box = $2

setenv IMAGIC_BATCH 1
echo "! "
echo "! "
echo "! ====================== "
echo "! IMAGIC ACCUMULATE FILE "
echo "! ====================== "
echo "! "
echo "! "
echo "! IMAGIC program: testim -----------------------------------------------"
echo "! "
/opt/qb3/imagic-110326/stand/testim.e <<EOF
junk,1,$tot
$box,$box
REAL
BLOBS
EOF
