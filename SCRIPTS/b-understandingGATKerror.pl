#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR inputs.pl
use strict;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# THE FILE PATH
my $filepath = $ARGV[0];
my $output = $ARGV[1];

#OPENING FILE
open (TABLE, "<$filepath") or die "Can't open file $filepath\n";
$/ = "##### ERROR   re"; 
my @fastaFILE = <TABLE>;
close(TABLE);

open (OUT, ">$output") or die "Can't open file $output\n";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
shift (@fastaFILE);
my $gencount = 0;
my $entry1 = $fastaFILE[0]; my @pieces1 = split(",", $entry1);
my $entry2 = $fastaFILE[1]; my @pieces2 = split(",", $entry2);
foreach my $tisk (0..$#pieces1){
	print OUT "$pieces1[$tisk]\t$pieces2[$tisk]\t$gencount\n";
	$gencount++;
}
$/ = "\n";  #returning back to the default
print "\n$gencount\n\n";
close(OUT);
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;

