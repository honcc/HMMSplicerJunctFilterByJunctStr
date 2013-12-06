#!/usr/bin/perl/ -w
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use List::Util qw ( sum );
use Storable;
use Data::Dumper::Names;

######################################################################################################################################################
#
#	Description
#		This is a perl script to extract the junctions in a HMMSplicer RAW BED file (i.e. 1 read per BED line) based on the a list of junction strings (i.e. cntg:start:end). It reads th BED file junctions and junctions strings 
#	, then it prints the junctions across multiple samples.
#
#	Input
#		--BEDPathListFile=		file path; path of the a file contains the list of the input BED files;
#		--junctStrListPath=		file path; path of a file contains a list of junctStr, one junctStr per line;
#		--outDir=				dir path: path of the output dir;
#
#	Usage
#		
#		perl HMMSplicerJunctFilterByJunctStr_v0.1.pl --BEDPathListFile=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/splicing/HMMSplicerJunctFilterByJunctStr/v0.1/HM1RahmanHMMSplicerBEDPaths.txt --junctStrListPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/splicing/HMMSplicerJunctFilterByJunctStr/v0.1/nonStochasticJunctStrNARPaperTableS3.txt
#
#	Assumption
#
#	History:
#		
#		v0.1
#		-debut
#
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#
#----------Read parameters ----------#
my ($BEDPathListFile, $junctStrListPath, $outDir) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

my ($junctstrHsh_ref) = getJunctStrList($junctStrListPath);

my ($BEDInfoHsh_ref) = readBEDFilePath($BEDPathListFile);
my $BEDDataHsh_ref;
($BEDDataHsh_ref, $BEDInfoHsh_ref) = readAllBEDFile($BEDInfoHsh_ref);
my %BEDInfoHsh = %{$BEDInfoHsh_ref};

filterBEDJunction($BEDDataHsh_ref, $junctstrHsh_ref, \%BEDInfoHsh, $outDir);

printCMDLogOrFinishMessage("finishMessage");

exit;
#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	my ($BEDPathListFile, $junctStrListPath, $outDir);
	
	$outDir = "./HMMSplicerJunctFilterByJunctStr/";

	foreach my $param (@ARGV) {
		if ($param =~ m/--BEDPathListFile=/) {$BEDPathListFile = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--junctStrListPath=/) {$junctStrListPath = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);} 
	}
	
	#---check file
	foreach my $fileToCheck ($BEDPathListFile, $junctStrListPath) {
		die "Cant read $fileToCheck" if not -s $fileToCheck;
		my @filePathAry = split /\//, $fileToCheck;
		my $fileName = $filePathAry[-1];
		print "$fileName checked.\n";
	}

	my @paramAry = ($BEDPathListFile, $junctStrListPath, $outDir);

	system "mkdir -p -m 777 $outDir/";
	open (PARAM, ">$outDir/parameters.txt");
	print PARAM Dumper($BEDPathListFile, $junctStrListPath, $outDir);
	close PARAM;
	
	return @paramAry;
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
########################################################################## readBEDFilePath
sub readBEDFilePath {

	#---my $BEDInfoHsh_ref = readBEDFilePath($BEDPathListFile);
	my $BEDPathListFile = $_[0];
	
	my %BEDInfoHsh;
	open (BEDPATHFILE, "$BEDPathListFile");
	while (my $theLine = <BEDPATHFILE>) {
		chomp $theLine;
		my ($sampleName, $replicate, $BEDPath) = split /\t/, $theLine;
		${${$BEDInfoHsh{$sampleName}}{$replicate}}{'BEDPath'} = $BEDPath;
		die "Cant read $BEDPath" if not -s $BEDPath;
		print "$sampleName $replicate BED checked.\n";
	}
	close BEDPATHFILE;
	
	return \%BEDInfoHsh;
}
########################################################################## readIndividualBEDFile
sub readIndividualBEDFile {

	#---my $BEDInfoHsh_ref = readIndividualBEDFile($BEDPath);
	
	my $BEDPath = $_[0];
	
	my %indivBEDDataHsh;
	my $totalReadNum = 0;
	open (BED, $BEDPath);
	while (my $theLine = <BED>) {
		chomp $theLine;
		next if ($theLine =~ m/^track name/);
		$totalReadNum++;
		my @theLineSplt = split (/\t/, $theLine);
		#DS571145	13916	14145	JUNC00000008	141	-	13916	14145	255,0,0	2	77,97	0,132
		my $cntg = $theLineSplt[0];
		my $bedStart = $theLineSplt[1];
		my $bedEnd = $theLineSplt[2];
		my @blkSizesSplt = split /,/, $theLineSplt[10];
		my $blk1Size = $blkSizesSplt[0];
		my $blk2Size = $blkSizesSplt[1];
		my $intronStart = $bedStart + $blk1Size + 1;
		my $intronEnd = $bedEnd - $blk2Size;
		my $junctStr = $cntg.":".$intronStart.":".$intronEnd; #---assumed to be unique
		${$indivBEDDataHsh{$junctStr}}{'readNum'}++;
	}
	close BED;
	
	my $junctNum = keys %indivBEDDataHsh;

	print "Totally $junctNum junctions read.\n";
	
	return \%indivBEDDataHsh, $totalReadNum;
}
########################################################################## readAllBEDFile
sub readAllBEDFile {
	
	#---my $BEDDataHsh_ref = readAllBEDFile($BEDInfoHsh_ref);

	my %BEDInfoHsh = %{$_[0]};
	
	my %BEDDataHsh;
	my @allTotalReadNumAry;
	foreach my $sampleName (sort {$a cmp $b} keys %BEDInfoHsh) {
		foreach my $replicate (sort {$a <=> $b} keys %{$BEDInfoHsh{$sampleName}}) {
			print "Reading BED of $sampleName $replicate.\n";
			my $BEDPath = ${${$BEDInfoHsh{$sampleName}}{$replicate}}{'BEDPath'};
			my ($indivBEDDataHsh_ref, $totalReadNum) = readIndividualBEDFile($BEDPath);
			${$BEDInfoHsh{$sampleName}{$replicate}}{'totalReadNum'} = $totalReadNum;
			push @allTotalReadNumAry, $totalReadNum;
			%{${$BEDDataHsh{$sampleName}}{$replicate}} = %{$indivBEDDataHsh_ref};
		}
	}
	
	my $avgTotalRdNum = sum(@allTotalReadNumAry)/@allTotalReadNumAry;
	
	foreach my $sampleName (sort {$a cmp $b} keys %BEDInfoHsh) {
		foreach my $replicate (sort {$a <=> $b} keys %{$BEDInfoHsh{$sampleName}}) {
			my $totalReadNum = ${$BEDInfoHsh{$sampleName}{$replicate}}{'totalReadNum'};
			my $scaleFactor = $avgTotalRdNum/$totalReadNum;
			${$BEDInfoHsh{$sampleName}{$replicate}}{'scaleFactor'} = $scaleFactor;
			print "$sampleName $replicate: totalReadNum = $totalReadNum, scaleFactor = $scaleFactor\n"
		}
	}

	return \%BEDDataHsh, \%BEDInfoHsh;
}
########################################################################## getGTFJunctStr
sub getJunctStrList {
	
	#---my ($GTFJunctHsh_ref, $GTFDataHsh_ref) = getGTFJunctStr($GTFPath);
	
	my $junctStrListPath = $_[0];
	
	my %junctListHsh;
	open (LIST, "$junctStrListPath");
	while (my $theLine = <LIST>) {
		chomp $theLine;
		$junctListHsh{$theLine}++;
	}
	
	my $junctNum = keys %junctListHsh;
	
	print "Totally $junctNum junctions stored from the list.\n";
	
	return \%junctListHsh;
}
########################################################################## filterBEDJunction
sub filterBEDJunction {
	
	my %BEDDataHsh = %{$_[0]};
	my %junctstrHsh = %{$_[1]};
	my %BEDInfoHsh = %{$_[2]};
	my $outDir = $_[3];
	
	my $blk1Size = my $blk2Size = 5;
	my %junctCountHsh;
	
	open (OUTPUT, ">$outDir/interSampleCount.log.xls");
	foreach my $junctStr (sort {$a cmp $b} keys %junctstrHsh) {
		my @outputAry;
		push @outputAry, "junctStr";
		foreach my $sampleName (sort {$a cmp $b} keys %BEDDataHsh) {
			foreach my $replicate (sort {$a <=> $b} keys %{$BEDDataHsh{$sampleName}}) {
				my $scaleFactor = ${$BEDInfoHsh{$sampleName}{$replicate}}{'scaleFactor'};
				$scaleFactor = sprintf "%.2f", $scaleFactor;
				push @outputAry, "$sampleName.$replicate [$scaleFactor]";
			}
		}
		push @outputAry, "totalReadNum";
		push @outputAry, "scaleReadNumSum";
		push @outputAry, "sampleAtLeast2Reads";
		print OUTPUT join "", ((join "\t", @outputAry), "\n");
		last;
	}

	foreach my $junctStr (sort {$a cmp $b} keys %junctstrHsh) {
		my @outputAry;
		push @outputAry, $junctStr;
		my $totalReadNum = 0;
		my $scaleReadNumSum = 0;
		my $sampleAtLeast2Reads = 0;
		foreach my $sampleName (sort {$a cmp $b} keys %BEDDataHsh) {
			foreach my $replicate (sort {$a <=> $b} keys %{$BEDDataHsh{$sampleName}}) {
				my $scaleFactor = ${$BEDInfoHsh{$sampleName}{$replicate}}{'scaleFactor'};
				my $readNum = 0;
				$readNum = ${${${$BEDDataHsh{$sampleName}}{$replicate}}{$junctStr}}{'readNum'} if exists ${${${$BEDDataHsh{$sampleName}}{$replicate}}{$junctStr}}{'readNum'};
				my $scaledReadNum = sprintf "%.2f", ($readNum*$scaleFactor);
				push @outputAry, $scaledReadNum;
				$totalReadNum += $readNum;
				$scaleReadNumSum += $scaledReadNum;
				$sampleAtLeast2Reads ++ if $readNum >= 2;
			}
		}
		push @outputAry, $totalReadNum;
		push @outputAry, $scaleReadNumSum;
		push @outputAry, $sampleAtLeast2Reads;
		print OUTPUT join "", ((join "\t", @outputAry), "\n");
	}
	close OUTPUT;
}
