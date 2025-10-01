#!/bin/bash
export CFOUR_TMPDIR=$(pwd)/tmp_cfour
mkdir -p $CFOUR_TMPDIR
CFOUR_TMPDIR=$CFOUR_TMPDIR xcfour > my_output_file.out
rm -rf $CFOUR_TMPDIR