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
my $fastq_dir = "";
my $work_dir = "";
my $bam_dir = "";
my $focus_dir = "";

my $qc_script = "run_qc.sh";
my $mapping_script = "run_mapping.sh";

if (@ARGV == 5) {
	$directory_file = $ARGV[0];
	$fastq_dir = $ARGV[1];
	$work_dir = $ARGV[2];
	$bam_dir = $ARGV[3];
	$focus_dir = $ARGV[4];
} else {
	print "Usage: perl configure_dapseq_qc_mapping.pl <directory file> <fastq directory> <working directory> <BAM output directory> <FOCUS output directory>\n";
	print "This script creates run_qc.sh and run_mapping.sh shell scripts for preprocessing and mapping of sequence reads for DAP-seq.\n";
	print "Scripts will be created in the working directory.\n";
	exit(0);
};


#Programs
my $extract_command = "zcat ";
my $trim_command = "SolexaQA++ dynamictrim ";
my $lengthsort_command = "SolexaQA++ lengthsort ";
my $align_command = "bowtie";
my $focus_command = "focus -q ";
my $samtools_view = "samtools view -bS -o";
my $samtools_sort = "samtools sort ";

#QC parameters
my $min_length = 30;

#Phylogenetic profiling parameters
my $focus_cutoff = "0.1";

#Mapping parameters
my $bowtie_libs_list = catfile($datadir,"libs.tsv");
my $m_option = "1";

#Peak calling parameters
my $q_value = "0.01";

###############

#Data structures
my %libs = ();
read_libs($bowtie_libs_list);

my %samples_directory = ();
read_directory($directory_file);


#my @samples = ();
#my @control_samples = ();

unless (-e $fastq_dir) {
	print "Directory $fastq_dir does not exist!\n";
	exit(1);
};

unless (-e $work_dir) {
	mkdir $work_dir;
};

unless (-e $bam_dir) {
	mkdir $bam_dir;
};

unless (-e $focus_dir) {
	mkdir $focus_dir;
};

$qc_script = catfile($work_dir, $qc_script);
$mapping_script = catfile($work_dir, $mapping_script);

open (QCOUTFILE, ">$qc_script") or die ("Unable to open file $qc_script");
print QCOUTFILE "#! /bin/sh\n\n";

open (MAPOUTFILE, ">$mapping_script") or die ("Unable to open file $mapping_script");
print MAPOUTFILE "#! /bin/sh\n\n";

opendir (DIR,"$fastq_dir");

my %fastq_files = ();
while (defined (my $file = readdir(DIR))) {
	if ($file =~ /R1_001\.fastq\.gz$/) { #don't forget to change
		$fastq_files{$file} = 1;
	}
}

foreach my $file (sort keys %fastq_files){
	process_single_fastq($file);
}

close DIR;
close QCOUTFILE;
close MAPOUTFILE;
chmod(0775, $qc_script);
chmod(0775, $mapping_script);

exit(0);

#######################
###   SUBROUTINES   ###
#######################

sub read_libs{
	my ($libs_file) = @_;
	open (INFILE, $libs_file) or die ("Unable to open libs file $libs_file");
	while (my $line = <INFILE>) {
		chomp $line;
		my ($lib_name, $lib_path, $size) = split(/\t/, $line);
		$libs{$lib_name}{"path"} = catfile($datadir,$lib_path);
		$libs{$lib_name}{"size"} = $size;
	}
}

sub read_directory{
	my $genome_name = "";
	my ($directory_file) = @_;
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
}

sub process_single_fastq {
	my ($file) = @_;
	my @filename_tokens = split(/_/, $file);
	my $sample = $filename_tokens[0]."_".$filename_tokens[1];
	unless (defined $samples_directory{$sample}){
		print "Sample for file $file not found in the directory.\n";
		return;
	}
	print QCOUTFILE "# Processing sample $sample\n";
	print QCOUTFILE $extract_command . $fastq_dir . "/". $file . " >" . $work_dir . "/" . $sample . ".fastq\n";
	print QCOUTFILE $trim_command . $work_dir . "/" . $sample . ".fastq\n";
	print QCOUTFILE "rm ". $work_dir . "/" . $sample . ".fastq\n";
	print QCOUTFILE $lengthsort_command . $work_dir . "/" . $sample . ".fastq.trimmed -l ". $min_length . "\n";
	print QCOUTFILE "rm ". $work_dir . "/" . $sample . ".fastq.trimmed\n";
	print QCOUTFILE "gzip ". $work_dir . "/" . $sample . ".fastq.trimmed.discard\n";
	print QCOUTFILE "mv ". $work_dir . "/" . $sample . ".fastq.trimmed.discard.gz ". $fastq_dir . "\n";
	print QCOUTFILE $focus_command . $work_dir . "/" . $sample . ".fastq.trimmed.single -m " . $focus_cutoff . "\n";
	print QCOUTFILE "gzip " . $work_dir . "/" . $sample . ".fastq.trimmed.single\n";
	print QCOUTFILE "mv " . $work_dir . "/" . $sample . ".fastq.trimmed.single__output.txt " . $focus_dir . "\n#######################################################\n\n";
	
	print MAPOUTFILE "# Processing sample $sample\n";
	print MAPOUTFILE "(echo $sample) >> mapping_quality.txt\n";
	print MAPOUTFILE "gunzip " . $work_dir . "/" . $sample . ".fastq.trimmed.single.gz\n";
	print MAPOUTFILE $align_command . " -S ". $libs{$samples_directory{$sample}{"genome"}}{"path"} . " -m " . $m_option . " " . $work_dir . "/" . $sample . ".fastq.trimmed.single " . $work_dir . "/" . $sample . ".sam 1>> mapping_quality.txt 2>&1\n";
	print MAPOUTFILE "gzip ". $work_dir . "/" . $sample . ".fastq.trimmed.single\n";
	print MAPOUTFILE "mv ". $work_dir . "/" . $sample . ".fastq.trimmed.single.gz ". $fastq_dir . "\n";
	print MAPOUTFILE "samtools view -bS -o " . $work_dir . "/" . $sample . ".bam " . $work_dir . "/" . $sample . ".sam\n";
	print MAPOUTFILE "rm " . $work_dir . "/" . $sample . ".sam\n";
	print MAPOUTFILE "samtools sort " . $work_dir . "/" . $sample . ".bam " . $bam_dir . "/" . $sample . ".sorted\n";
	print MAPOUTFILE "rm " . $work_dir . "/" . $sample . ".bam\n#######################################################\n\n";

}
