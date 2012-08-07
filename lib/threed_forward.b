#!/bin/csh -f

set file = $1
set ang = $2

setenv IMAGIC_BATCH 1
echo "! "
echo "! "
echo "! ====================== "
echo "! IMAGIC ACCUMULATE FILE "
echo "! ====================== "
echo "! "
echo "! "
echo "! IMAGIC program: forward_3d -------------------------------------------"
echo "! "
/opt/qb3/imagic-110326/threed/forward_3d.e PROJECTION_MODE FORWARD <<EOF
$file
my_forw
-99999
YES
ASYM_TRIANGLE
C1
EQUIDIST
ZERO
$ang
NO
EOF
/opt/qb3/imagic-110326/incore/incnorvar.e <<EOF
SOF
my_forw
my_forw1
0.9,0.2
10.0
EOF

rm my_forw.*
mv my_forw1.img my_forw.img
mv my_forw1.hed my_forw.hed
rm my_forw1.*
