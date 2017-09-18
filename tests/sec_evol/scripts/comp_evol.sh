#!/bin/bash

usage () {
        echo " Usage: pre_evol.sh <astname>"
#        echo " comments ..."
        exit 1
}
if [ $# -lt 1 ]; then
        echo 'no astname supplied'
        usage
fi

rm -f comp_evol_fil.m;
rm -f comp_evol_dat.m;
rm -f sec2fla.pl; 
rm -f dat2fla.x;
rm -f fil2fla.x;
ln -s ../mlab/comp_evol_fil.m .;
ln -s ../mlab/comp_evol_dat.m .;
ln -s ../scripts/sec2fla.pl .;
ln -s ../dat2fla.x .;
ln -s ../fil2fla.x .;

sec2fla.pl $1;
dat2fla.x;
fil2fla.x;
