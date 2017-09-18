#!/usr/bin/perl -w

$ast=$ARGV[0];
$clones=$ARGV[1];
#chomp $clones;

for ($i=1;$i<=$clones;$i++){
#    print $i;
    open(INFILE,"<clones.opt");
    open(OUTFILE,">orb9.opt.$ast\_$i");

    $count=0;
    while($riga=<INFILE>){
	if($count==8){
	    print OUTFILE "$ast\_$i\n";
	}else{
	    print OUTFILE $riga;
	}
	$count++;
    }

    close(OUTFILE);
    close(INFILE);

}
