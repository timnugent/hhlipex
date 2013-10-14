#!/usr/bin/perl -w
# Predict transmembrane lipid exposure

use strict;
use Getopt::Long;
use File::Basename;

## Set these two paths:
my $hhpath = "/home/tnugent/Applications/hhsuite-2.0.13-linux-x86_64/bin";
my $dbname = "/home/tnugent/databases/uniprot20_02Sep11";

my $svm_model = "model/lipex.model";
my $output_path = "output";

sub usage{

	print "Usage:\n";
	print "$0 [-p <hhblits bin path>] [-d <hhblit database] [-o <output path>] <fasta file>\n";
	exit;
}

sub get_arguments{

	# Process command line arguments
	if (!$ARGV[0]){
		&usage;
	}else{
		my $result = GetOptions ("d=s" => \$dbname,
					"p=s" => \$hhpath,
					"m=s" => \$svm_model,
					"o=s" => \$output_path,
					"h"  => sub {&usage;});
	}
}

sub main{

	&get_arguments();

	my $fa = $ARGV[-1];
	die "$fa doesn't exist!\n" unless -e $fa;
	my @suffix = qw(fa fasta);
	my($filename,$directories,$suffix) = fileparse($fa,qr/\.\D.*/);
	my $aln = $output_path."/".$filename.".aln";
	my $hhm = $output_path."/".$filename.".hhm";
	my $out = $output_path."/".$filename.".out";
	my $lipex = $output_path."/".$filename.".lipex";

	my $hhblits = "$hhpath/hhblits";
	die "$hhblits doesn't exist!\n" unless -e $hhblits;
	my $hhmake  = "$hhpath/hhmake";
	die "$hhmake doesn't exist!\n" unless -e $hhmake;
	my $hhpssm = "bin/hhpssm";
	die "$hhpssm doesn't exist! Run make\n" unless -e $hhpssm;

	my $command = "$hhblits -pcm 2 -pca 0.9 -i $fa -oa3m $aln -n 3 -e 1e-6 -d $dbname";	
	print "$command\n" unless -e $aln;
	`$command` unless -e $aln;
	$command = "$hhmake -pcm 2 -v 4 -i $aln -o $hhm > $out";	
	print "$command\n";
	`$command`;
	$command = "$hhpssm -m $svm_model $out > $lipex";
	print "$command\n";
	print `$command`;
	print "\nTransmembrane lipid exposure probability estimates have been written to:\n$lipex\n\n";

}
&main;

