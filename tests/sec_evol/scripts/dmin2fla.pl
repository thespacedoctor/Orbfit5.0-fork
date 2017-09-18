#!/usr/bin/perl -w

$ast=$ARGV[0];
$secdat=$ARGV[1]; #1=sec; 2=dat

open(INF,"<$ast.dmin");
open(OUTF,">dmin.fla");
$count=0;
while($riga=<INF>){
    $count++;
    if($count==1){
	@pezzi=split(/\s+/,$riga);
	$tast0=$pezzi[2];
#	print OUTF "$tast0  0  0  0  0  0\n";
    }elsif($count>2){
	print OUTF "$riga";
    }
}

close(OUTF);
close(INF);

if($secdat==1){
    system "echo $tast0  0  0  0  0  0 > dmin.fla.tmp";
}elsif($secdat==2){
    system "echo $tast0  0  0  0  > dmin.fla.tmp";
}

system "cat dmin.fla >> dmin.fla.tmp";
system "cp -p dmin.fla.tmp dmin.fla";
