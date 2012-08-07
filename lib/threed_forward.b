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
