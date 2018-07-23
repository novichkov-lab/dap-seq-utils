#! /usr/bin/perl
use 5.010;
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use File::Spec::Functions 'catfile';

#Paths and file names
my $datadir = catfile(dirname(dirname(abs_path($0))),"data");
my $directory_file = "";
my $work_dir = "";
my $bam_dir = "";
my $output_dir = "";

if (@ARGV == 4) {
	$directory_file = $ARGV[0];
	$work_dir = $ARGV[1];
	$bam_dir = $ARGV[2];
	$output_dir = $ARGV[3];
} else {
	print "Usage: perl configure_dapseq_peak_calling.pl <directory file> <working directory> <BAM files directory> <output directory>\n";
	print "This script creates run_macs_meme.sh shell script for peak calling and motif search.\n";
	print "Scripts will be created in the working directory.\n";
	exit(0);
};


my $genome_library_file = catfile($datadir,"libs.tsv");
my $gbk_list = catfile($datadir,"genome_list.txt");
my $macs_script = "run_macs_meme.sh";
my $log_file = catfile($work_dir, "log.txt");


# MACS2 parameters
my $genome_size = "";
my $extsize = "250";
my $qvalue = "0.0001";
my $keepdup = "auto";
my $noextsize = 1;

# Parameters for peak filtering before MEME search
my $enrichment_cutoff = 3;
my $pvalue_cutoff = 10;
my $qvalue_cutoff = 10;

###############
my %libs = ();
my %samples_directory = ();


#my @samples = ();
#my @control_samples = ();

unless (-e $work_dir) {
	mkdir $work_dir;
};
open (LOGFILE, ">$log_file") or die ("Unable to open file $log_file");

# read genomes library
read_libs($genome_library_file);

# read directory file
read_directory($directory_file);

unless (-e $bam_dir) {
	print LOGFILE "Directory $bam_dir does not exist!\n";
	exit(1);
}
unless (-e $output_dir) {
	mkdir $output_dir;
};

$macs_script = catfile($work_dir, $macs_script);
open (OUTFILE, ">$macs_script") or die ("Unable to open file $macs_script");
print OUTFILE "#! /bin/sh\n\n";

opendir (DIR,"$bam_dir");
while (defined (my $file = readdir(DIR))) {
	if ($file =~ /\.sorted\.bam$/) {
		process_bam_file($file);
	}
}
close DIR;

close OUTFILE;
chmod(0775, $macs_script);
print "$macs_script created \n";
close LOGFILE;
exit(0);

###################
### SUBROUTINES ###
###################

sub read_libs{
	my ($libs_file) = @_;
	open (INFILE, $libs_file) or die ("Unable to open libs file $libs_file");
	while (my $line = <INFILE>) {
		chomp $line;
		my ($lib_name, $lib_path, $size) = split(/\t/, $line);
		$libs{$lib_name}{"path"} = catfile($datadir,$lib_path);
		$libs{$lib_name}{"size"} = $size;
	}
	return 1;
}

sub read_directory{
	my ($directory_file) = @_;
	my $genome_name = "";
	open (INFILE, $directory_file) or die ("Unable to open libs file $directory_file");
	while (my $line = <INFILE>) {
		chomp $line;
		if ($line =~ /^#/) {
			next;
		}
		my ($sample, $genome, $protein, $replicate, $treatment, $control) = split(/\t/, $line);
		unless ($genome_name){
			$genome_name = $genome;
		}
		unless ($genome eq $genome_name){
			print "Different genome names in the directory file\n";
			die();
		}
		$samples_directory{$sample}{"genome"} = $genome;
		$samples_directory{$sample}{"protein"} = $protein;
		$samples_directory{$sample}{"replicate"} = $replicate;
		$samples_directory{$sample}{"treatment"} = $treatment;
		$samples_directory{$sample}{"control"} = $control;
	}
	return 1;
}

sub process_bam_file {
	my ($file) = @_;
	my $sample = $file;
	$sample =~ s/\.sorted\.bam$//;
	unless (defined $samples_directory{$sample}){
		print LOGFILE "Sample for file $file not found in the directory.\n";
		return 0;
	}

	unless (exists $libs{$samples_directory{$sample}{"genome"}}){
		print LOGFILE "Genome " . $samples_directory{$sample}{"genome"} . "not found in the library\n";
		return 0;
	}

	my $control_bam_file = catfile($bam_dir, $samples_directory{$sample}{"control"} . ".sorted.bam");
	unless (-e $control_bam_file){
		print LOGFILE "Control file $control_bam_file for sample $sample not found.\n";
		return 0;
	}
	print OUTFILE "\n#######################################################\n";
	print OUTFILE "macs2 callpeak -t " . $bam_dir . "/" . $file . " -c " . $control_bam_file . " -f BAM -g " . $libs{$samples_directory{$sample}{"genome"}}{"size"} . " -n " . $output_dir . "/" . $sample . "_vs_" . $samples_directory{$sample}{"control"} . " -q " . $qvalue;
	unless ($noextsize){
		print OUTFILE " --extsize " . $extsize;
	}
	print OUTFILE " --nomodel --keep-dup=" . $keepdup . "\n";
#	print OUTFILE "perl /mnt/data2/DAP-seq/scripts/run_meme_on_peaks_v3.pl " . $output_dir . "/" . $sample . "_vs_" . $samples_directory{$sample}{"control"} . "_peaks.xls " . $output_dir . " \"" . $samples_directory{$sample}{"genome"} . "\" -1\n";
	print OUTFILE "perl /mnt/data2/DAP-seq/scripts/run_meme_on_filtered_peaks.pl " . $output_dir . "/" . $sample . "_vs_" . $samples_directory{$sample}{"control"} . "_peaks.xls " . $output_dir . " \"" . $samples_directory{$sample}{"genome"} . "\" " . $enrichment_cutoff . " " . $pvalue_cutoff . " " . $qvalue_cutoff . "\n";
	print OUTFILE "#######################################################\n";
	return 1;
}
