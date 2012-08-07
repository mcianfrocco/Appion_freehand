#!/bin/csh -f

set file = $1
set align = $2
set out = $3
set micro = $4
set def = $5
set mag = $6

setenv IMAGIC_BATCH 1
echo "! "
echo "! "
echo "! ====================== "
echo "! IMAGIC ACCUMULATE FILE "
echo "! ====================== "
echo "! "
echo "! "
echo "! IMAGIC program: imagic_frealign --------------------------------------"
echo "! "
/opt/qb3/imagic-110326/frealign/imagic_frealign.e MODE CONVERT <<EOF
$file
$align
$out
MRA
I
$mag
PER
$micro
PER-IM
$def
CTF
EOF
