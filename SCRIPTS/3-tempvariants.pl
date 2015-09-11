#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR ...

=pod

=head1 NAME

NONE 

=head1 SYNOPSIS

xxx.pl [--help] [--manual]

=head1 DESCRIPTION

Accepts all folders from frnakenstein output.
 
=head1 OPTIONS

=over 3

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   DBI
   DBD::mysql
   Getopt::Long
   Pod::Usage

=head1 AUTHOR

Written by Modupe Adetunji, 
Center for Bioinformatics and Computational Biology Core Facility, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2015 MOA.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's usage
=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR 
use strict;
use File::Basename;
use DBI;
use DBD::mysql;
use Getopt::Long;
use Pod::Usage;

#CREATING LOG FILES
my $std_out = '/home/modupeore17/.LOG/allivar-'.`date +%m-%d-%y_%T`; chomp $std_out; $std_out = $std_out.'.log';
my $std_err = '/home/modupeore17/.LOG/allivar-'.`date +%m-%d-%y_%T`; chomp $std_err; $std_err = $std_err.'.err';
my $jobid = "ALLI-".`date +%m-%d-%y_%T`;
my $progressnote = "/home/modupeore17/.LOG/allinote".`date +%m-%d-%y_%T`; chomp $progressnote; $progressnote = $progressnote.'.txt'; 

open(STDOUT, '>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$std_err") or die "Error file doesn't exist";
 
#ARGUMENTS
my($help,$manual,$in1);
GetOptions (	
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);

#file path for input THIS SHOULD BE CONSTANT
$in1 = "/opt/apache2/frankdec/subdirectories/mapcount_output";
#making sure the input file is parsable
my @temp = split('',$in1); $in1 = undef; my $checking = pop(@temp); push (@temp, $checking); unless($checking eq "/"){ push (@temp,"/")}; foreach(@temp){$in1 .= $_};

# DATABASE ATTRIBUTES
my $dsn = 'dbi:mysql:transcriptatlas';
my $user = 'frnakenstein';
my $passwd = 'maryshelley';

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

my $mystderror = "Contact Modupe Adetunji amodupe\@udel.edu\n";

# FOLDER VARIABLES
my ($top_folder, $cuff_folder, $htseq_folder);

# RESULTS_HASH
my %Hashresults; my %Birdresults; my %HTSEQ;

# DATABASE VARIABLES
my ($dbh, $sth, $syntax, $row, @row);
my $code="DE"; 

#DIRECTORY
my (@parse, @NewDirectory);

# TABLE VARIABLES
my ($parsedinput, $len);
my ($lib_id, $total, $mapped, $unmapped, $deletions, $insertions, $junctions, $genes, $isoforms, $prep, $date); #transcripts_summary TABLE
my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat); # GENES_FPKM & ISOFORMS_FPKM TABLE

#VARIANTS FOLDER & HASHES
my $Mfolder; my %VCFhash; my %VEPhash; 

#PARSABLE GENOMES FOR ANALYSIS
my $GENOMES="/home/modupeore17/.GENOMES/";
my %parsablegenomes = ("chicken" => 1,"alligator" => 2, ); #genomes that work.


#INDEPENDENT PROGRAMS
my $PICARDDIR="/home/modupeore17/.software/picard-tools-1.136/picard.jar";
my $GATKDIR="/home/modupeore17/.software/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar";
my $VEP="/home/modupeore17/.software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl";
 

#GETTING ACCURATE FILE PATHS : recommended but not needed
system `locate $in1 > parse.txt`;
my $filelocation = `head -n1 parse.txt`; @parse = split('\/', $filelocation); $in1 = undef; $len =$#parse-1;foreach(0..$len){$in1 .= $parse[$_]."\/";};
system `rm -f parse.txt`;

#OPENING FOLDER
opendir(DIR,$in1) or die "Folder \"$in1\" doesn't exist\n"; 
my @Directory = readdir(DIR);
close(DIR);
#pushing each subfolder
foreach (@Directory){
	if ($_ !~ /^\w*_\d*$/ && $_ =~ /^[a-zA-Z].*$/){
		push (@NewDirectory, $_);
	}
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#CREATING EMAIL NOTIFICATION
my $notification = '/home/modupeore17/.LOG/note.txt';
my $email = 'amodupe@udel.edu';
open (NOTE, ">$notification");
print NOTE "Subject: Starting job : $jobid\n\nName of log files\n\t$std_out\n\t$std_err\n";
system "sendmail $email < $notification";
close NOTE;
system "rm -rf $notification";
# -----------------------------------

# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE : $dsn\n\n";
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";
#DELETENOTDONE();
foreach my $NewFolder (@NewDirectory) { 
	my $input = $in1.$NewFolder;
	opendir(NEW,$input) or die "Folder \"$input\" doesn't exist\n";
	my @Subdirectory = readdir(NEW);
	close(NEW);
	foreach my $SubNewFolder (@Subdirectory) {
		if ($SubNewFolder =~ /^\w.*_(\d.*)$/){
			CHECKING();
			if($1 == 1052){#unless (exists $Hashresults{$1}){
				if (exists $Birdresults{$1}){
					$parsedinput = "$input/$SubNewFolder";
		            $Mfolder = "$parsedinput/variant_output";
					PARSING($1,$parsedinput,$NewFolder);
					#progress report
					open (NOTE, ">>$progressnote");
					print NOTE "Subject: Update notes : $jobid\n\nCompleted library\t$1\n";
					system "sendmail $email < $progressnote"; close NOTE;
				}
				else {
					print "\nSkipping \"library_$1\" in \"$NewFolder\" folder because it isn't in birdbase\n$mystderror\n";
				}
			}
		}	 
	}
}
SUMMARYstmts();
system "rm -rf $progressnote";
# DISCONNECT FROM THE DATABASE
print "\n\tDISCONNECTING FROM THE DATABASE : $dsn\n\n";
$dbh->disconnect();

#send finish notification
my $endnotification = '/home/modupeore17/.LOG/endnote.txt';
open (NOTE, ">$endnotification");
print NOTE "Subject: Job completed : $jobid\n\nName of log files\n\t$std_out\n\t$std_err\n";
system "sendmail $email < $endnotification";
close NOTE;
system "rm -rf $endnotification";
# -----------------------------------
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - -S U B R O U T I N E S- - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub DELETENOTDONE {
	print "\n\tDELETING NOT DONE\n\n";
	#CHECKING TO MAKE SURE NOT "done" FILES ARE REMOVED
	$syntax = "select library_id from transcripts_summary where status is NULL";
	$sth = $dbh->prepare($syntax);
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	my $incompletes = undef; my $count=0; my @columntoremove;
	while ($row = $sth->fetchrow_array() ) {
		$count++;
		$incompletes .= substr($row,0,-2).", ";
	}
	if ($count >= 1){
		$incompletes = substr($incompletes,0,-2);         
		#DELETE FROM variants_result
		$syntax = "delete from variants_result where library_id in \( $incompletes \)";
		$sth = $dbh->prepare($syntax); $sth->execute();
                #DELETE FROM variants_summary
                $sth = $dbh->prepare("delete from variants_summary where library_id in ( $incompletes )"); $sth->execute();
		#DELETE FROM genes_fpkm
		$sth = $dbh->prepare("delete from genes_fpkm where library_id in ( $incompletes )"); $sth->execute();
		#DELETE FROM isoforms_fpkm
		$sth = $dbh->prepare("delete from isoforms_fpkm where library_id in ( $incompletes )"); $sth->execute();
		#DELETE FROM frnak_metadata
		$sth = $dbh->prepare("delete from frnak_metadata where library_id in ( $incompletes )"); $sth->execute();
		#DELETE FROM transcripts_summary
		$sth = $dbh->prepare("delete from transcripts_summary where library_id in ( $incompletes )"); $sth->execute();
	}
}
sub CHECKING {
	#CHECKING THE LIBRARIES ALREADY IN THE DATABASE
	$syntax = "select library_id from transcripts_summary";
	$sth = $dbh->prepare($syntax);
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	my $number = 0;
	while ($row = $sth->fetchrow_array() ) {
    	my $newrow = substr($row,0,-2);
		$Hashresults{$newrow} = $number; $number++;
	}
	$syntax = "select library_id from bird_libraries";
	$sth = $dbh->prepare($syntax);
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	$number = 0;
	while ($row = $sth->fetchrow_array() ) {
    	my $newrow = substr($row,0,-2);
		$Birdresults{$newrow} = $number; $number++;
	}
}
sub PARSING {
	print "\n\tINSERTING TRANSCRIPTS INTO THE DATABASE : \t library_$_[0]\n\n";
	my $first = $_[0]; my $second = $_[1];
	$lib_id = "$first$code";
	$top_folder = "$second/tophat_out"; $cuff_folder = "$second/cufflinks_out"; $htseq_folder = "$second/htseq_output";
	
	#created a log file check because
	#I'm having issues with the authenticity of the results caz some of the log files is different
	my $logfile = `head -n 1 $top_folder/logs/run.log`;
	my @allgeninfo = split('\s',$logfile);

	#also getting metadata info
	my ($str, $ann, $ref, $seq) = 0; #defining calling variables
	#making sure I'm working on only the chicken files for now, need to find annotation of alligator
	my @checkerno = split('\/',$allgeninfo[$#allgeninfo]);
	my @numberz =split('_', $checkerno[$#checkerno]);
	if ($numberz[0] == $first){
		#making sure the arguments are accurately parsed
		if ($allgeninfo[1] =~ /.*library-type$/){$str = 2; $ann = 4; $ref = 9; $seq = 10;}
		elsif($allgeninfo[3] =~ /\-o$/){$str=99; $ann=99; $ref = 5; $seq = 6;}
		else {print "File format doesn't match what was encoded for frnak_metadata\nSyntax error:\n\t@allgeninfo\n$mystderror\n";next;}
		#assuring we are working on available genomes
		my $refgenome = (split('\/', $allgeninfo[$ref]))[-1]; #reference genome name
		if (exists $parsablegenomes{$refgenome}){
			open(ALIGN, "<$top_folder/align_summary.txt") or die "Can't open file $top_folder/align_summary.txt\n";
 			open(GENES, "<$cuff_folder/genes.fpkm_tracking") or die "Can't open file $cuff_folder/genes.fpkm_tracking\n";
 			open(ISOFORMS, "<$cuff_folder/isoforms.fpkm_tracking") or die "Can't open file $cuff_folder/isoforms.fpkm_tracking\n";
 		
 	# PARSER FOR transcripts_summary TABLE
			while (<ALIGN>){
				chomp;
				if (/Input/){my $line = $_; $line =~ /Input.*:\s+(\d+)$/;$total = $1;}
				if (/Mapped/){my $line = $_; $line =~ /Mapped.*:\s+(\d+).*$/;$mapped = $1;}
			} close ALIGN;
			$unmapped = $total-$mapped;
			$deletions = `cat $top_folder/deletions.bed | wc -l`; $deletions--;
			$insertions = `cat $top_folder/insertions.bed | wc -l`; $insertions--;
			$junctions = `cat $top_folder/junctions.bed | wc -l`; $junctions--;
			$genes = `cat $cuff_folder/genes.fpkm_tracking | wc -l`; $genes--;
			$isoforms = `cat $cuff_folder/isoforms.fpkm_tracking | wc -l`; $isoforms--;
			$prep = `cat $top_folder/prep_reads.info`;
			$date = `date +%Y-%m-%d`;

		#PARSING FOR SNPanalysis
			my $accepted = "$top_folder/accepted_hits.bam"; @parse = split('\/\/',$accepted); $accepted = undef; $len = $#parse+1; foreach(@parse){$accepted .= $_; if($len>1){$accepted .="\/"; $len--;}};
			my $run_log = "$top_folder/logs/run.log"; my @parse = split('\/\/',$run_log); $run_log = undef; $len = $#parse+1; foreach(@parse){$run_log .= $_; if($len>1){$run_log .="\/"; $len--;}};
		
			#INSERT INTO DATABASE : transcriptatlas
			#transcripts_summary table
			# $sth = $dbh->prepare("insert into transcripts_summary (library_id, total_reads, mapped_reads, unmapped_reads, deletions, insertions, junctions, isoforms, genes, info_prep_reads, date ) values (?,?,?,?,?,?,?,?,?,?,?)");
# 			$sth ->execute($lib_id, $total, $mapped, $unmapped, $deletions, $insertions, $junctions, $isoforms, $genes, $prep, $date);
			#GENES_FPKM table
# 			$sth = $dbh->prepare("insert into genes_fpkm (library_id, tracking_id, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id, chrom_no, chrom_start, chrom_stop, length, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
		# 	while (<GENES>){
# 				chomp;
# 				my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
# 				unless ($track eq "tracking_id"){ #check & specifying undefined variables to null
# 					if($class =~ /-/){$class = undef;} if ($ref_id =~ /-/){$ref_id = undef;}
# 					if ($length =~ /-/){$length = undef;} if($coverage =~ /-/){$coverage = undef;}
# 					$locus =~ /^(.+)\:(.+)\-(.+)$/;
# 					$chrom_no = $1; $chrom_start = $2; $chrom_stop = $3;
# 					$sth ->execute($lib_id, $track, $class, $ref_id, $gene, $gene_name, $tss, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat );
# 				}
# 			} close GENES;

			#ISOFORMS_FPKM table
# 			$sth = $dbh->prepare("insert into isoforms_fpkm (library_id, tracking_id, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id, chrom_no, chrom_start, chrom_stop, length, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
			# while (<ISOFORMS>){
# 				chomp;
# 				my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
# 				unless ($track eq "tracking_id"){
# 					if ($class =~ /-/){$class = undef;} if ($ref_id =~ /-/){$ref_id = undef;}
# 					if ($length =~ /-/){$length = undef;} if($coverage =~ /-/){$coverage = undef;}
# 					$locus =~ /^(.+)\:(.+)\-(.+)$/;
# 					$chrom_no = $1; $chrom_start = $2; $chrom_stop = $3;
# 					$sth ->execute($lib_id, $track, $class, $ref_id, $gene, $gene_name, $tss, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat );
# 				}
# 			} close ISOFORMS;
# 		
			#extracting the syntax used for METADATA 
			#temporary replace the file path
			my ($replace, $temptest,$stranded, $sequences, $annotationfile, $annfileversion);
			unless ($ann ==99){
				$temptest = `head -n 1 $allgeninfo[$ann]`;
				if (length($temptest) < 1 ){
					$replace = "/home/modupeore17/.big_ten"; foreach my $change (0..$#allgeninfo){if ($allgeninfo[$change] =~ /^(.*big_ten).*$/){$allgeninfo[$change] =~ s/$1/$replace/g;}}
				}
				$annotationfile = uc ( (split('\.',((split("\/", $allgeninfo[$ann]))[-1])))[-1] ); #(annotation file)
                                $annfileversion = substr(`head -n 1 $allgeninfo[$ann]`,2,-1); #annotation file version
        		}
			else { $annotationfile = "NULL"; $annfileversion = "NULL"; }
			if ($str == 99){ $stranded = "NULL"; }else { $stranded = $allgeninfo[$str]; } # (stranded or not)	
			my $otherseq = $seq++;
			unless(length($allgeninfo[$otherseq])<=1){ #sequences 
				$sequences = ( ( split('\/', $allgeninfo[$seq]) ) [-1]).",". ( ( split('\/', $allgeninfo[$otherseq]) ) [-1]);
			} else { $sequences = ( ( split('\/', $allgeninfo[$seq]) ) [-1]);}
			
			#frnak_metadata table
# 			$sth = $dbh->prepare("insert into frnak_metadata (library_id,ref_genome, ann_file, ann_file_ver, stranded, sequences,user ) values (?,?,?,?,?,?,?)");
# 			$sth ->execute($lib_id, $refgenome, $annotationfile, $annfileversion, $stranded,$sequences,$_[2] );
			#print "$lib_id, $refgenome, $annotationfile, $annfileversion, $stranded,$sequences,$_[2]\n";
			#htseq table: Has to be created on the fly for different genomes; NOT WORKING (08/20/2015) too many columns
			#HTSEQ($htseq_folder,$first,$refgenome);
		
			#variant analysis
			VARIANTS($first, $accepted, $refgenome);

			#Finally : the last update. transcripts_summary table updating status column with 'done'
# 			$sth = $dbh->prepare("update transcripts_summary set status='done' where library_id = $first");
# 			$sth ->execute();
		}
		else { 
			my $parsabletemp = 1;
            print "The reference genome isn't available, available genomes are : ";
            foreach my $pargenomes (keys %parsablegenomes){
                print "\"$pargenomes\"";
                if($parsabletemp < (keys %parsablegenomes)){ print ", "; $parsabletemp++; }
            }
            print " rather what you have is $allgeninfo[$ref]\n$mystderror\n";
		}
	}
	else {print "library_id dont match $numberz[0] == $first\n";}
}
sub HTSEQ {
	print "\n\tSTORING HTSEQ IN THE DATABASE\n\n";

# too many columns, I have to figure out another way to import this table: maybe NoSQL.
	my %HTSEQ = undef; my $syntaxforhtseq = undef; my %dbtables = undef;
	open(HTSEQ, "<$_[0]/$_[1]".".counts") or die "Can't open file $_[0]/$_[1]".".counts\n";
	while (<HTSEQ>){
		chomp;
		my ($NAME, $VALUE) = split /\t/;
		$HTSEQ{uc($NAME)} = $VALUE;
	} close HTSEQ;
	
	#---------------------
	#show tables in the databases;
	$syntax = "show tables"; my $dbcount=0;
	$sth = $dbh->prepare($syntax);
	$sth->execute or die "SQL Error : (me) SHOWING TABLES : $DBI::errstr\n";
	while ($row = $sth->fetchrow_array() ) {
		$dbtables{$row} = $dbcount++;
	}
	my $prospectivetable = $_[2]."_HTSEQ";
	unless (exists $dbtables{$prospectivetable}){my $i = 0;
		my $startofsyntax = "create table $prospectivetable (library_id varchar(50) not null, ";
		foreach (sort keys %HTSEQ){$i++;
			if (/^[0-9a-zA-Z]/){
				$syntaxforhtseq .= "`$_` int, ";
			}
		}			
		my $endofsyntax = "primary key (library_id), foreign key (library_id) references transcripts_summary (library_id));";
		$syntax = $startofsyntax.$syntaxforhtseq.$endofsyntax; print "number of columns $i++\n\n";
		$sth = $dbh->prepare($syntax);
		$sth->execute or die "SQL Error : (me) Creating HTSEQ table: $DBI::errstr\n";
	}
	#------------------------------
	my ($namesyntax, $valuesyntax, $questionsyntax) = undef; my $keyscount = 0;
	foreach (keys %HTSEQ){
		$keyscount++;
		$namesyntax .= "$_, ";
		$valuesyntax .= "$HTSEQ{$_}, ";
		$questionsyntax .= "?,";
	}
	$namesyntax = substr($namesyntax,0,-2);
	$valuesyntax = substr($valuesyntax,0,-2);
	$questionsyntax = substr($questionsyntax,0,-1);
	$sth = $dbh->prepare("insert into $prospectivetable ( $namesyntax ) values ($questionsyntax)");
	$sth-> execute($valuesyntax);
}
sub VARIANTS{
	print "\n\tWORKING ON VARIANT ANALYSIS\n\n";
 	my $libraryNO = "library_".$_[0]; my $bamfile = $_[1]; 
	my $REF= "/home/modupeore17/.big_ten/alligator/alligator/alligator.fa";
	my $specie = $_[2];
 	my $DICT = "$GENOMES/$_[2]/$_[2]".".dict";

	#Making directory & testing the existence
	`mkdir $Mfolder`;
	my $testpath = `ls $Mfolder`;
	if (length($testpath) < 1 ){
		$Mfolder = "/home/modupeore17/CHICKENvariants/$libraryNO"; `mkdir $Mfolder`;
	}

	##PICARD
#	`java -jar $PICARDDIR SortSam INPUT=$bamfile OUTPUT=$Mfolder/$libraryNO.bam SO=coordinate`;
		
	#ADDREADGROUPS
#	my $addreadgroup = "java -jar $PICARDDIR AddOrReplaceReadGroups INPUT=$Mfolder/$libraryNO.bam OUTPUT=$Mfolder/$libraryNO"."_add.bam SO=coordinate RGID=LAbel RGLB=Label RGPL=illumina RGPU=Label RGSM=Label";
#	`$addreadgroup`;
	
	#MARKDUPLICATES
#	my $markduplicates = "java -jar $PICARDDIR MarkDuplicates INPUT=$Mfolder/".$libraryNO."_add.bam OUTPUT=$Mfolder/".$libraryNO."_mdup.bam M=$Mfolder/".$libraryNO."_mdup.metrics CREATE_INDEX=true";
#	`$markduplicates`;
	
	#SPLIT&TRIM
	my $splittrim = "java -jar $GATKDIR -T SplitNCigarReads -R $REF -I $Mfolder/".$libraryNO."_mdup.bam -o $Mfolder/".$libraryNO."_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar";
	`$splittrim`;
	
	#GATK
	my $gatk = "java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $Mfolder/".$libraryNO."_split.bam -o $Mfolder/$libraryNO.vcf";
	`$gatk`;
	
	#perl to select DP > 5 & get header information
	FILTERING($Mfolder, "$Mfolder/$libraryNO.vcf");

	#ANNOTATIONS : running VEP
	my $veptxt = "perl $VEP -i $Mfolder/".$libraryNO."_DP5.vcf --fork 24 --species $specie  --cache --merged --everything on --terms ensembl --output_file $Mfolder/".$libraryNO."_VEP.txt";
	`$veptxt`;
	my $vepvcf = "perl $VEP -i $Mfolder/".$libraryNO."_DP5.vcf --fork 24 --species $specie  --cache --vcf --merged --everything on --terms ensembl --output_file $Mfolder/".$libraryNO."_VEP.vcf";
	`$vepvcf`;
	
	#DATABASE INSERT
# 	DBVARIANTS($Mfolder."/".$libraryNO."_VEP.vcf", $libraryNO);
} 
sub FILTERING {
	my $input = $_[1];
	my $wkdir = $_[0];
	unless(open(FILE,$input)){
		print "File \'$input\' doesn't exist\n";
		exit;
	}

	my $out = fileparse($input, qr/(\.vcf)?$/);

	my $output = "$out"."_DP5.vcf";
	open(OUT,">$wkdir/$output");
	my $output2 = "$out"."_header.vcf";
	open(OUT2,">$wkdir/$output2");

	my @file = <FILE>; chomp @file; close (FILE);
	
	foreach my $chr (@file){
		if ($chr =~ /^chr/){
			my @chrdetails = split('\t', $chr);
			my $chrIwant = $chrdetails[7];
			my @morechrsplit = split(';', $chrIwant);
			foreach my $Imptchr (@morechrsplit){
				if ($Imptchr =~ m/^DP/) {
					my @addchrsplit = split('=', $Imptchr);
					if ($addchrsplit[1] > 4){print OUT "$chr\n";}
				}
			}
		}
		else {
			print OUT "$chr\n"; print OUT2 "$chr\n";
		}
	}
	close (OUT); close (OUT2);
}
sub DBVARIANTS{	
	print "\n\tINSERTING VARIANTS INTO THE DATABASE\n\n";
	#disconnecting and connecting again to database just incase
	$dbh->disconnect(); $dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";

	$_[1] =~ /^library_(\d*)$/;
	my $libnumber = "$1$code";
	my $folder = undef;

	#VEP file
	my @splitinput = split('\/', $_[0]);
	foreach my $i (0..$#splitinput-1){$folder.="$splitinput[$i]/";$i++;}
	my $information = fileparse($_[0], qr/(\.vcf)?$/);
	
	#variant_metadata
	 my $gatk_version = ( ( split('\/',$GATKDIR)) [-2] ); my $vep_version = ( ( split('\/',$VEP)) [-4] ); my $picard_version = ( ( split('\/',$PICARDDIR)) [-2] );

	#OUTPUT file
	my $output = "$folder$information".".table";
	open(OUTDBVAR,">$output");
	print OUTDBVAR "#CHR\tPOS\tREF\tALT\tQUAL\tCLASS\t";
	print OUTDBVAR "ENSEMBL-GENE\tENTREZ-GENE\tSYMBOL\tTRANSCRIPT\t";
	print OUTDBVAR "FEATURE\tTYPE OF GENE\tCONQ\tPROTEIN\t";
	print OUTDBVAR "PROTEIN POSITION\tCODON CHANGE\tZYGOSITY\tdbSNP\n";
	my $date = `date +%Y-%m-%d`;

#initializing the hash tables . . .
	%VCFhash = ();
	%VEPhash = ();
#running through subroutines . . . 
	VEPVARIANT($_[0]);
	my ($itsnp,$itindel,$itvariants) = 0;
#printing to output table & variant_results table
#VARIANT_SUMMARY
	$sth = $dbh->prepare("insert into variants_summary ( library_id, VEP_version, Picard_version, GATK_version, date ) values (?,?,?,?,?)");
	$sth ->execute($libnumber, $vep_version, $picard_version, $gatk_version, $date);
#VARIANT_RESULTS
	foreach my $abc (sort keys %VCFhash) {
		foreach my $def (sort {$a <=> $b} keys %{ $VCFhash{$abc} }) {
			my @vcf = split('\|', $VCFhash{$abc}{$def});
			my @vep = split('\|', $VEPhash{$abc}{$def});
			if ($vcf[3] =~ /,/){
				my $first = split(",",$vcf[3]);
				if (length $vcf[2] == length $first){
					$itvariants++; $itsnp++;
				}
				else {
					$itvariants++; $itindel++;
				}
			}
			elsif (length $vcf[2] == length $vcf[3]){
				$itvariants++; $itsnp++;
			}
			else {
				$itvariants++; $itindel++;
			}
			print OUTDBVAR "$abc\t$def\t$vcf[2]\t$vcf[3]\t$vcf[4]\t$vep[0]\t";
			print OUTDBVAR "$vep[2]\t$vep[11]\t$vep[10]\t$vep[3]\t";
			print OUTDBVAR "$vep[1]\t$vep[9]\t$vep[4]\t$vep[5]\t";
			print OUTDBVAR "$vep[6]\t$vep[7]\t$vcf[5]\t$vep[8]\n";
		

			$sth = $dbh->prepare("insert into variants_result ( library_id, chrom, position, ref_allele, alt_allele, quality, variant_class,
			ensembl_gene, entrez_gene, gene_name, transcript, feature, gene_type,
			consequence, protein, protein_position, codon_change, zygosity, existing_variant ) values
			(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
			$sth ->execute($libnumber, $abc, $def, $vcf[2], $vcf[3], $vcf[4], $vep[0], $vep[2], $vep[11], $vep[10], $vep[3], $vep[1], $vep[9], 
			$vep[4], $vep[5], $vep[6], $vep[7], $vcf[5], $vep[8]);
		}
	}
	close (OUTDBVAR);
	#VARIANT_SUMMARY
	$syntax = "update variants_summary set total_VARIANTS = $itvariants, total_SNPS = $itsnp, 
	total_INDELS = $itindel, status = \'done\' where library_id like \"$libnumber\"";
	$sth = $dbh->prepare($syntax);
	$sth ->execute();
}
sub VEPVARIANT {
	#working on VEP variants
	my %Geneinfo = ''; my %Transcriptinfo = ''; my %Conqinfo = ''; my %Varclass = '';
	my %Featureinfo = ''; my %Proinfo = ''; my %Prochangeinfo = ''; my %Codoninfo = '';
	my %dbSNPinfo = ''; my %locinfo = ''; my %Entrezinfo = ''; my %GENEtype = '';
	my %GENEname = ''; my %GeneEntrezinfo = ''; my $position; my %TranscriptEntrezinfo = '';
	my %location = '';
	unless(open(FILE,$_[0])){print "File \'$_[0]\' doesn't exist\n";exit;}
	my $verd;
	my @file = <FILE>;
	chomp @file;
	close (FILE);
	foreach my $chr (@file){
		unless ($chr =~ /^#/){
			my @chrdetails = split('\t', $chr);
			
			#removing the random chromosomes (scaffolds) - because no vital information can be found for them.
			my @morechrsplit = split(';', $chrdetails[7]);
			if (((split(':', $chrdetails[9]))[0]) eq '0/1'){$verd = "heterozygous";}
			elsif (((split(':', $chrdetails[9]))[0]) eq '1/1'){$verd = "homozygous";}
			elsif (((split(':', $chrdetails[9]))[0]) eq '1/2'){$verd = "heterozygous alternate";}
	
			#VCFhash information
			$VCFhash{$chrdetails[0]}{$chrdetails[1]} = "$chrdetails[0]|$chrdetails[1]|$chrdetails[3]|$chrdetails[4]|$chrdetails[5]|$verd";
			
			#Processing the VEP section
			my @finalchrsplit = split("\,",(((split('=',((split(';',$chrdetails[7]))[-1])))[1]))); 
			foreach my $FCR (0..$#finalchrsplit){
				my @vepdetails = split('\|', $finalchrsplit[$FCR]);	
				$location{"$chrdetails[0]|$chrdetails[1]"} = "$chrdetails[0]|$chrdetails[1]"; #specifying location.
				if ($vepdetails[1] !~ /WITHIN_NON_CODING_GENE/){
					if ($vepdetails[4] =~/^ENSG.*/){
						#GENE - 1
						if (exists $Geneinfo{$chrdetails[0]}{$chrdetails[1]}){
							if ($Geneinfo{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[4]){
								$Geneinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[4];
							}
							else {
								my $tempdetails = "$Geneinfo{$chrdetails[0]}{$chrdetails[1]},$vepdetails[4]";
								$Geneinfo{$chrdetails[0]}{$chrdetails[1]} = $tempdetails;
							}
						}
						else {	
							$Geneinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[4];
						}
					}
					else {
						#GENE - 1 - ENTREZ
						if (exists $GeneEntrezinfo{$chrdetails[0]}{$chrdetails[1]}){
							if ($GeneEntrezinfo{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[4]){
								$GeneEntrezinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[4];
							}
							else {
								my $temp0details = "$GeneEntrezinfo{$chrdetails[0]}{$chrdetails[1]},$vepdetails[4]";
								$GeneEntrezinfo{$chrdetails[0]}{$chrdetails[1]} = $temp0details;
							}
						}
						else {	
							$GeneEntrezinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[4];
						}
					}
					#TRANSCRIPT - 2
					if ($vepdetails[6] =~/^ENSG.*/){
						if (exists $Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}){
							if ($Transcriptinfo{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[6]){
								$Transcriptinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[6];
							}
							else {
								my $temp2details = "$Transcriptinfo{$chrdetails[0]}{$chrdetails[1]},$vepdetails[6]";
								$Transcriptinfo{$chrdetails[0]}{$chrdetails[1]} = $temp2details;
							}
						}
						else {	
							$Transcriptinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[6];
						}
					}
					else { #ENTREZ
						if (exists $TranscriptEntrezinfo{$chrdetails[0]}{$chrdetails[1]}){
							if ($TranscriptEntrezinfo{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[6]){
								$TranscriptEntrezinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[6];
							}
							else {
								my $temp2bdetails = "$TranscriptEntrezinfo{$chrdetails[0]}{$chrdetails[1]},$vepdetails[6]";
								$TranscriptEntrezinfo{$chrdetails[0]}{$chrdetails[1]} = $temp2bdetails;
							}
						}
						else {	
							$TranscriptEntrezinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[6];
						}
					}
				
					#CONSEQUENCE - 3
					if (exists $Conqinfo{$chrdetails[0]}{$chrdetails[1]} ){
						if ($Conqinfo{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[1]){
							$Conqinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[1];
						}
						else {
							my $temp3details = "$Conqinfo{$chrdetails[0]}{$chrdetails[1]},$vepdetails[1]";
							$Conqinfo{$chrdetails[0]}{$chrdetails[1]} = $temp3details;
						}			
					}
					else {	
						$Conqinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[1];
					}
							
					#VARIANT CLASS - 4
					if (exists $Varclass{$chrdetails[0]}{$chrdetails[1]}){
						if ($Varclass{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[20]){
						$Varclass{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[20];
						}
						else {
							my $temp4details = "$Varclass{$chrdetails[0]}{$chrdetails[1]},$vepdetails[20]";
							$Varclass{$chrdetails[0]}{$chrdetails[1]} = $temp4details;
						}			
					}
					else {	
						$Varclass{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[20];
					}		
						
					#FEATURE TYPE - 5
					if (exists $Featureinfo{$chrdetails[0]}{$chrdetails[1]}){
						if ($Featureinfo{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[5]){
							$Featureinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[5];
						}
						else {
							my $temp5details = "$Featureinfo{$chrdetails[0]}{$chrdetails[1]},$vepdetails[5]";
							$Featureinfo{$chrdetails[0]}{$chrdetails[1]} = $temp5details;
						}			
					}
					else {	
						$Featureinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[5];
					}
					
					#PROTEIN POSITION - 6
					if (exists $Proinfo{$chrdetails[0]}{$chrdetails[1]}){
						if ($Proinfo{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[14]){
							$Proinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[14];
						}
						else {
							my $temp6details = "$Proinfo{$chrdetails[0]}{$chrdetails[1]},$vepdetails[14]";
							$Proinfo{$chrdetails[0]}{$chrdetails[1]} = $temp6details;
						}			
					}
					else {	
						$Proinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[14];
					}
		
					#PROTEIN CHANGE - 7 (or Amino Acid)
					if (exists $Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}){
						if ($Prochangeinfo{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[15]){
							$Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}= $vepdetails[15];
						}
						else {
							my $temp7details = "$Prochangeinfo{$chrdetails[0]}{$chrdetails[1]},$vepdetails[15]";
							$Prochangeinfo{$chrdetails[0]}{$chrdetails[1]} = $temp7details;
						}			
					}
					else {	
						$Prochangeinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[15];
					}
					
					#CODON CHANGE - 8
					if (exists $Codoninfo{$chrdetails[0]}{$chrdetails[1]}){
						if ($Codoninfo{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[16]){
							$Codoninfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[16];
						}
						else {
							my $temp8details = "$Codoninfo{$chrdetails[0]}{$chrdetails[1]},$vepdetails[16]";
							$Codoninfo{$chrdetails[0]}{$chrdetails[1]} = $temp8details;
						}			
					}
					else {	
						$Codoninfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[16];
					}
		
					#dbSNP - 9
					if (exists $dbSNPinfo{$chrdetails[0]}{$chrdetails[1]}){
						if ($dbSNPinfo{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[17]){
							$dbSNPinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[17];
						}
						else {
							my $temp9details = "$dbSNPinfo{$chrdetails[0]}{$chrdetails[1]},$vepdetails[17];";
							$dbSNPinfo{$chrdetails[0]}{$chrdetails[1]} = $temp9details;
						}			
					}
					else {	
						$dbSNPinfo{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[17];
					}			
					
					#GENE name - 10
					if (exists $GENEname{$chrdetails[0]}{$chrdetails[1]}){
						if ($GENEname{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[3]){
							$GENEname{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[3];
						}
						else {
							my $temp10details = "$GENEname{$chrdetails[0]}{$chrdetails[1]},$vepdetails[3]";
							$GENEname{$chrdetails[0]}{$chrdetails[1]} = $temp10details;
						}
					}
					else {	
						$GENEname{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[3];
					}
					
					#GENE type - 11
					if (exists $GENEtype{$chrdetails[0]}{$chrdetails[1]}){
						if ($GENEtype{$chrdetails[0]}{$chrdetails[1]} eq $vepdetails[7]){
							$GENEtype{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[7];
						}
						else {
							my $temp11details = "$GENEtype{$chrdetails[0]}{$chrdetails[1]},$vepdetails[7]";
							$GENEtype{$chrdetails[0]}{$chrdetails[1]} = $temp11details;
						}
					}
					else {	
						$GENEtype{$chrdetails[0]}{$chrdetails[1]} = $vepdetails[7];
					}
				}
			}
		}
	}

	foreach my $alldetails (keys %location){
		my ($chrdetails1, $chrdetails2) = split('\|', $alldetails);
		#cleaning up the text
		my $clean1 = CLEANUP($Varclass{$chrdetails1}{$chrdetails2});
		my $clean2 = CLEANUP($Featureinfo{$chrdetails1}{$chrdetails2});
		my $clean3 = CLEANUP($Geneinfo{$chrdetails1}{$chrdetails2});
		my $clean4 = CLEANUP($Transcriptinfo{$chrdetails1}{$chrdetails2});
		my $clean5 = CLEANUP($Conqinfo{$chrdetails1}{$chrdetails2});
		my $clean6 = CLEANUP($Proinfo{$chrdetails1}{$chrdetails2});
		my $clean7 = CLEANUP($Prochangeinfo{$chrdetails1}{$chrdetails2});
		my $clean8 = CLEANUP($Codoninfo{$chrdetails1}{$chrdetails2});
		my $clean9 = CLEANUP($dbSNPinfo{$chrdetails1}{$chrdetails2});
		my $clean10 = CLEANUP($GENEtype{$chrdetails1}{$chrdetails2});
		my $clean11 = CLEANUP($GENEname{$chrdetails1}{$chrdetails2});
		my $clean12 = CLEANUP($GeneEntrezinfo{$chrdetails1}{$chrdetails2});
		
		$VEPhash{$chrdetails1}{$chrdetails2} = "$clean1|$clean2|$clean3|$clean4|$clean5|$clean6|$clean7|$clean8|$clean9|$clean10|$clean11|$clean12";
	}
}
sub CLEANUP {
	#cleaning up the VEP variants so that it doesn't have repetitions in the output
	my @unclean = split(',', $_[0]);
	my ($cleansyntax, %Hashdetails, %Hash2) = undef;
	foreach (0..$#unclean){
		if ($unclean[$_] =~ /^[a-zA-Z0-9]/){
			$Hashdetails{$unclean[$_]} = $_;
		}
	}
	foreach my $unique (keys %Hashdetails){
		if ($unique =~ /^[a-zA-Z0-9]/){
			$Hash2{$Hashdetails{$unique}} = $unique;
		}
	}
	foreach my $final (sort keys %Hash2){
		if ($Hash2{$final} =~ /^[a-zA-Z0-9]/){
			$cleansyntax .= ",$Hash2{$final}";
		}
	}
	my $returnclean = substr($cleansyntax, 1);
	return $returnclean;
}
sub SUMMARYstmts{
	#transcripts_summary from the database
	print "\n\tEXECUTING SELECT STATEMENT ON THE DATABASE TABLES \n";
	$syntax = "select count(*) from transcripts_summary";
	$sth = $dbh->prepare($syntax);
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"transcripts_summary\" table \t:\t @row\n";}
	#genes_fpkm
	$syntax = "select count(*) from genes_fpkm";
	$sth = $dbh->prepare($syntax);
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"genes_fpkm\" table \t\t:\t @row\n";}
	#isoforms_fpkm
	$syntax = "select count(*) from isoforms_fpkm";
	$sth = $dbh->prepare($syntax);
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"isoforms_fpkm\" table \t:\t @row\n";}
	#variant_summary
	$syntax = "select count(*) from variants_summary";
	$sth = $dbh->prepare($syntax);
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"variants_summary\" table \t:\t @row\n";}
	#variant_list
	$syntax = "select count(*) from variants_result";
	$sth = $dbh->prepare($syntax);
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"variants_results\" table \t:\t @row\n";}
	#frnak_metadata
	$syntax = "select count(*) from frnak_metadata";
        $sth = $dbh->prepare($syntax);
        $sth->execute or die "SQL Error: $DBI::errstr\n";
        while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"frnak_metadata\" table \t:\t @row\n";}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
close STDOUT; close STDERR;
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;
