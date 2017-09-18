#!/bin/bash

usage () {
        echo " Usage: VAsecevol.sh <astname>"
#        echo " comments ..."
        exit 1
}
if [ $# -lt 1 ]; then
        echo 'no astname supplied'
        usage
fi

cp -p ./VAs $1/;
#cp -p ../lostPHAs/$1/*.eq1 epoch/
rm -f epoch/$1\_*.eq1;
cp -p $1\_*.eq1 epoch/;

rm -f $1/$1\_*.eq1;
mv -f $1\_*.eq1 $1/;

FILES=$1/$1\_*.eq1;

i=0;
for f in $FILES
do
    i=$(($i+1));
    echo "Processing $1_$i file"
#    cp -p $f epoch/$1\_$i.eq1; # questo comando crea casino
  
    echo $1_$i | sec_evol.x;

#    mv -f $1\_$i.sec $1/$1\_$i.sec;
    mv -f $1\_$i.pre $1/$1\_$i.pre;
    mv -f $1\_$i.pro $1/$1\_$i.pro;

done
echo "$1: processed $i VAs"

