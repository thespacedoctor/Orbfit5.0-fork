#!/usr/bin/perl -w

open(INF,"<vastfla");
open(OUF,">vastflanew");

$count=1;
while($riga=<INF>){
    $num=2*int($count/2);
#    print "$num    $count\n";
    if($num==$count){
#	print "$riga\n";
	print OUF "$riga";
    }
    $count++;
}

close(OUF);
close(INF);
