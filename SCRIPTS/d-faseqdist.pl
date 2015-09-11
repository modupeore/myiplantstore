#! usr/bin/perl


#OBJECTIVE

print "\n\tFrequency of chromosomal sequence length \n";

my $input = $ARGV[0];
my $output = $ARGV[1];
my $i=0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}
open (OUT, ">$output");
$/ = ">";
while (<FILE>){
	chomp $_;
	my @test = split("\n", $_);
	my $header = $test[0];
	my $seq = length($test[1]);
	print OUT "$header\t$seq\n";
}
print "\nDone\n";
close (FILE); close (OUT);

