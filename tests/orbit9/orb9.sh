#!/bin/bash

usage () {
        echo " Usage: orb9.sh <astname> <nclones>"
#        echo " comments ..."
        exit 1
}
if [ $# -lt 1 ]; then
        echo 'no astname supplied'
        usage
fi
if [ $# -lt 2 ]; then
        echo 'no nclones supplied'
        usage
fi

if [ $2 -lt 1 ]; then
    echo " no. clones should be > 0 "
    exit 1;
fi

mv -f orb9.opt orb9.opt.tmp;

# eliminate mercury and uranus from the initial input file 'allplamxx.inc'
# there is a problem with the time, it must be the same of inbaric.opt and 
# the same of asteroid initial epoch
inbarmerc.x;   #produce barmerxx.inc allplmxx.inc
pvensat.pl $3; #generate pvensat.inc from allplmxx.inc
cp -p pvensat.inc clones/pvensat.inc;

cd clones;

# create a catalog of clones for the asteroid
FILE=$1.cat
if [ ! -f $FILE ]; then
    echo "File $FILE not found!";
    exit 1;
fi
FILE=clon_cat.pl
if [ ! -f $FILE ]; then
    echo "File $FILE not found!";
    exit 1;
fi
rm -f $1.cat.clones;
echo "... clones/$1.cat.clones";
clon_cat.pl $1 $2;

# create a catalog of clones for the planets
FILE=pla_clones.pl
if [ ! -f $FILE ]; then
    echo "File $FILE not found!";
    exit 1;
fi
rm -f pvensat_*.inc;
echo "... clones/pvensat_<nclone>.inc";
pla_clones.pl $2;

# create an option file with the name of the clone
FILE=clones.opt
if [ ! -f $FILE ]; then
    echo " File $FILE not found ! ";
    exit 1;
fi
rm -f orb9.opt.$1_*;
echo "... clones/orb9.opt.<clone>";
orbopt.pl $1 $2;

cd ..;

# create a link to the asteroid clone catalog
rm -f ast.cat.clones;
ln -s clones/$1.cat.clones ast.cat.clones;

# remove previous outputs
rm -f vast_*;
rm -f vpla_*;

n=1;
for ((i=1;i<=$2;i++)); do
    echo "$1_$i";
# create a link to the asteroid clone option file
    rm -f orb9.opt;
    ln -s clones/orb9.opt.$1_$i orb9.opt;
    for ((j=1;j<=$2;j++)); do
# create a link to the planets clone option file
	rm -f pvensat.inc
	ln -s clones/pvensat_$j.inc pvensat.inc;
# *************** ORBIT9.x **************************
	echo "orbit9.x $n";
	time orbit9.x;
	mv vpla.dat vpla_$n.dat;
	mv vast.dat vast_$n.dat;
	mv vpla.fil vpla_$n.fil;
	mv vast.fil vast_$n.fil;
# **************************************************
	n=$[n+1];
    done;
done;

rm -f ast.cat.clones;
rm -f orb9.opt;
rm -f pvensat.inc;

#
cp -p orb9.opt.tmp orb9.opt;

echo " All done ! "
