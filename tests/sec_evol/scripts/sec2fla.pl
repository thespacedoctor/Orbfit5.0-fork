#!/usr/bin/perl -w

$ast=$ARGV[0];

open(INF,"<$ast.sec");
open(OUTF,">secfile.fla");
$count=0;
while($riga=<INF>){
    $count++;
    if($count==1){
	@pezzi=split(/\s+/,$riga);
	$tast0=$pezzi[1];
	$aa=$pezzi[2];
	print OUTF "$tast0  $aa  0  0  0  0  0  0\n";
    }else{
	print OUTF "$riga";
    }
}

close(OUTF);
close(INF);
