#!/bin/bash

usage () {
        echo " Usage: comp_circ.sh <astname>"
#        echo " comments ..."
        exit 1
}
if [ $# -lt 1 ]; then
        echo 'no astname supplied'
        usage
fi

rm -f comp_circ.m;
rm -f mev2fla.pl; 
ln -s ../mlab/comp_circ.m .;
ln -s ../scripts/mev2fla.pl .;

cp -f ../../propneo2/$1.mev $1.circ.mev;
mev2fla.pl $1.circ;
mv mevfile.fla mevcirc.fla

mev2fla.pl $1;

# Plot asteroid seculare evolution (mevfile)
# and asteroid coplanar circular secular evolution (mevcirc)
#MATLAB -> COMP_CIRC

#rm -f suplabel.m;
#ln -s ../mlab/suplabel.m suplabel.m;
#matlab -nodesktop -nosplash -r "run comp_circ;"
