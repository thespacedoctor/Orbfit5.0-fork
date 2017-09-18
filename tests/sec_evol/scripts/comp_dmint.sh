#!/bin/bash

usage () {
    echo " Usage: comp_moid.sh <astname>"
#   echo " comments ..."
    exit 1
}
if [ $# -lt 1 ]; then
    echo 'no astname supplied'
    usage
fi

(echo sec ;echo $1) |dmint_evol.x;
cp -p $1.dmin $1/dmint.sec;
dmin2fla.pl $1 1;
mv dmin.fla $1/dminsec;
(echo dat ;echo $1) |dmint_evol.x;
cp -p $1.dmin $1/dmint.dat;
dmin2fla.pl $1 2;
mv dmin.fla $1/dmindat;

rm -f $1.dmin;

#matlab -nodesktop -nosplash -r "run moid_433;" #quit;"

#gv moid.eps
#gv moid_bw.eps