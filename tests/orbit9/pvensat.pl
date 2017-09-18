#!/usr/bin/perl -w

$epoch=$ARGV[0];

open(INF,"<allplmxx.inc");
open(OUTF,">pvensat.inc");

print OUTF "   5 VENUS     EARTHMOON MARS      JUPITER   SATURN\n";
print OUTF "t0=   2.4",$epoch,"5000000000D+06 ;  JPL DE 405  23 LUG 2010   0.0000\n";
print OUTF "coosys='KEP HEL ECLM00 DEG';ilce=   0\n";
print OUTF "  a(AU)    e    I    Node  nrev   arg.peri nrev  an.med.  nrev\n";
print OUTF "; planetary inverse masses as used\n";
print OUTF "  1.0000000000000002D+00  4.0852371000000008D+05  3.2890056140000006D+05\n";
print OUTF "  3.0987080000000088D+06  1.0473485999999953D+03  3.4978979999996586D+03\n";
print OUTF "  0.00000000\n";


$count=0;
while($riga=<INF>){
    $count++;
#    print "COUNT=$count\n";
    if($count>=13 && $count<=17){
	print OUTF "$riga";
    }   
}

close(OUTF);
close(INF);
