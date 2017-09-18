#!/bin/bash

usage () {
        echo " Usage: dmin_evol.sh <astname>"
#        echo " comments ..."
        exit 1
}
if [ $# -lt 1 ]; then
        echo 'no astname supplied'
        usage
fi

FILES=$1_*.sec
#echo $FILES
#exit;

cd $1;
rm -f dmin2fla.pl;
ln -s ../scripts/dmin2fla.pl;
cd ..;

i=0;
for f in $FILES
do
#    echo $f;
    i=$(($i+1));    
    (echo sec; echo $1_$i) | dmint_evol.x;
    mv -f $1_$i.dmin $1;
    cd $1;
    dmin2fla.pl $1_$i 1;
    mv dmin.fla dmin_$i.fla;
    cd ..;
done
# cleaning
mv -f $1_*.sec $1/;
echo "$1: processed $i VAs"

