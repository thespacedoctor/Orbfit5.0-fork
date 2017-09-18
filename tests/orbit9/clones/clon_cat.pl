#!/usr/bin/perl -w


$ast=$ARGV[0];
$nclones=$ARGV[1];

open(INFILE,"<$ast.cat");
open(OUTFILE,">$ast.cat.clones");

$count=1;
while($riga=<INFILE>){
    if($count<=6){
	print OUTFILE $riga;
    }else{
	@campi=split(/\s+/,$riga);
	chomp @campi;
	$nclo=1;
	while($nclo<$nclones+1){
	    $duepi=360.0;
	    $campo7=($nclo-1)*$duepi/$nclones;
	    printf OUTFILE "%7s%8s%12s%3s%22s%3s%22s%3s%22s%3s%22s%3s%22s%3s%22.16e%1s%5s%2s%4s\n","\'".$ast._.$nclo."\'","",$campi[1],"",$campi[2],"",$campi[3],"",$campi[4],"",$campi[5],"",$campi[6],"",$campo7,"",$campi[8],"",$campi[9] ;
	$nclo++;
	}
    }
    $count++;
}

close(OUTFILE);
close(INFILE);
