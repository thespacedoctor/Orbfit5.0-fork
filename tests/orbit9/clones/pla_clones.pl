#!/usr/bin/perl -w

$n=$ARGV[0];
$duepi=360.0;

for ($i=1;$i<=$n;$i++){
#    print "qui \$i vale $i e \$n vale $n\n";
    open(INFILE,"<pvensat.inc");
    open(OUTFILE,">pvensat_$i.inc");

    $ang=($i-1)*$duepi/$n;
    $count=1;

    while($riga=<INFILE>){
	if($count<=8){
	    print OUTFILE $riga;
	}else{
	    @campi=split(/\s+/,$riga);
	    chomp @campi;
	    printf OUTFILE "%2s%12s%1s%12s%3s%10s%2s%11s%2s%1s%2s%11s%2s%1s%2s%12.8f%3s%1s\n","",$campi[1],"",$campi[2],"",$campi[3],"",$campi[4],"",$campi[5],"",$campi[6],"",$campi[7],"",$ang,"",$campi[9];
	}
	$count++;
    }
    close(OUTFILE);
    close(INFILE);
}


