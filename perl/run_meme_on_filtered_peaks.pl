#!/usr/bin/perl

use 5.010;
use strict;
use warnings;
use Cwd qw();
use File::Basename;
use File::Spec::Functions 'catfile';

my $cwd = Cwd::cwd();
my $datadir = dirname($cwd);
$datadir = catfile($datadir,"data");
my $infile = "";
my $output_dir = "";
my $genome_list_file = catfile($datadir,"genome_list.txt");

my $genome_name = "";
my $limit = -1; #obsolete
my $enrichment_cutoff = 0;
my $pvalue_cutoff = 0;
my $qvalue_cutoff = 0;

my $max_fragment_length = 12000;
my $exclude_peaks_within_genes = 0; 
#MEME parameters
my @meme_cmd=('meme');
my $mod = "anr";
my $nmotifs = "1";
my $minw = "12";
my $maxw = "28";
 
# Data structures
my %genomes_list = ();
my %qvalue_peaks_HoA = (); #hash key is -log10(qvalue), value is an array of peak data lines


my $header = "chr	start	end	length	abs_summit	pileup	-log10(pvalue)	fold_enrichment	-log10(qvalue)	name";

if (@ARGV == 6) {
	$infile = $ARGV[0];
	$output_dir = $ARGV[1];
	$genome_name = $ARGV[2];
	$enrichment_cutoff = $ARGV[3];
	$pvalue_cutoff = $ARGV[4];
	$qvalue_cutoff = $ARGV[5];
} else {
	print "Usage: perl run_meme_on_peaks.pl <input file name with list of MACS2 peaks without extension> <directory where output files will be written> <genome name> <enrichment cut-off> <-log10(pvalue) cutoff> <-log10(qvalue) cutoff>\n";
	print "This script creates fasta file with sequences of peaks found by MACS and runs MEME with this file\n";
	exit(0);
};

print "Parameters: enrichment cutoff " . $enrichment_cutoff . " pvalue cutoff " . $pvalue_cutoff . " qvalue cutoff " . $qvalue_cutoff . "\n";

#read least of sequences 
read_genome_list($genome_list_file);
unless (exists $genomes_list{$genome_name}){
	print "Genome $genome_name not found in the $genome_list_file file.\n";
	exit(1);
}

#read sequences
my %sequences = ();
my %genes = ();
my %upstreams = ();

foreach my $seqfile (@{$genomes_list{$genome_name}}){
	if (-e $seqfile){
		print $seqfile."\n";
		processGbkFile($seqfile) ;
	} else {
		print "ERROR: File $seqfile not found!\n";
		exit(1);
	}
}

if (($output_dir ne "")&&(!(-e $output_dir))) {
	print "Directory $output_dir does not exist!\n";
	exit(1);
} elsif ($output_dir eq "") {
	$output_dir = ".";
};

#$infile = $output_dir . "/" . $infile;

open (INFILE, "$infile") or die ("Cannot open file $infile");
while (my $line = <INFILE>) {
	chomp $line;
	if (($line !~ /^\t.*/) && ($line !~ /^#.*/) && ($line ne "") && ($line !~ /^\Q$header/)){
		my @entry = split (/\t/, $line);
		if (($entry[6] >= $pvalue_cutoff) && ($entry[7] >= $enrichment_cutoff) && ($entry[8] >= $qvalue_cutoff)) {
#			print "\nPeak found: " .$line;
			if (exists $qvalue_peaks_HoA{$entry[8]}){
				push @{$qvalue_peaks_HoA{$entry[8]}}, $line;
			} else {
				my @arr = ($line);
				$qvalue_peaks_HoA{$entry[8]} = \@arr;
			}
		}
	};
}
close INFILE;

my $outfile = catfile($output_dir,"peaks.fasta");
open (OUTFILE, ">$outfile") or die ("Cannot open file $outfile");
foreach my $qvalue (sort {$b <=> $a} keys %qvalue_peaks_HoA){
	if ($limit == 0){
		last;
	}
	foreach my $line (@{$qvalue_peaks_HoA{$qvalue}}){
		$line = get_peak_sequence ($line);
		if ($line) {
			$limit--;
			print OUTFILE $line."\n";		
		}
	}
}
close OUTFILE;

my $meme_outfile = $infile . "_motifs_$mod.txt";

push @meme_cmd, $outfile;
push @meme_cmd, '-o';
push @meme_cmd, $output_dir;
push @meme_cmd, '-text';
push @meme_cmd, '-dna';
push @meme_cmd, '-mod';
push @meme_cmd, $mod;
push @meme_cmd, '-nmotifs';
push @meme_cmd, $nmotifs;
push @meme_cmd, '-minw';
push @meme_cmd, $minw;
push @meme_cmd, '-maxw';
push @meme_cmd, $maxw;
push @meme_cmd, '-revcomp';
push @meme_cmd, '-pal';
push @meme_cmd, '-maxsize 1000000';
push @meme_cmd, '>' . $meme_outfile;

print join(" ", @meme_cmd);

system (join(" ", @meme_cmd));
exit(0);

#######################
##### SUBROUTINES #####
#######################


sub read_genome_list{
	my ($file) = @_;
	open (INFILE, $file) or die ("Unable to open libs file $file");
	while (my $line = <INFILE>) {
		chomp $line;
		unless ($line !~ /^#/){
			next;
		}
		unless ($line ne ""){
			next;
		}
		my (undef, $genome, $filepath) = split(/\t/, $line);
		if (exists $genomes_list{$genome}){
			push @{$genomes_list{$genome}}, $filepath;
		} else {
			my @arr = ($filepath);
			$genomes_list{$genome} = \@arr;
		}
	}
	return 1;
}

sub processGbkFile {
	my ($seqfile) = @_;
	my $seq_flag = 0;
	my $accession = "";
	my @sequence = ();
	my $gene_flag = 0;
	my $gene_feature = "^     CDS";
	my $gene_qualifier = "                     \/gene=\"";
	my $locus_tag_qualifier = "^                     \/locus_tag=\"";
	my $old_locus_tag_qualifier = "^                     \/old_locus_tag=\"";
	my $feature_start = 0;
	my $feature_end = 0;
	my $feature_strand = "not set";
	my $size = 0;
	my @genes_from_gbk = ();
	my @upstreams_data = ();


	
	open (SEQFILE, "$seqfile");
	while (my $line = <SEQFILE>) {
		chomp $line;
		if ($line =~ /^ORIGIN.*/) {
			if ($accession eq "") {
				print "ERROR: Accession line not found in file $seqfile \nCheck file format\n";
				close SEQFILE;
				exit(1);
			}
			$seq_flag = 1;
#		} elsif ($line =~ /^ACCESSION   /) {
#			$line =~ s/^ACCESSION   //g;
#			$accession = $line;
		} elsif ($line =~ /^VERSION     /) {
			$line =~ s/^VERSION     //g;
			($accession) = split(/\s/,$line);
		} elsif ($seq_flag) {
			push @sequence, process_seq_line($line);
		} elsif ($line =~ /$gene_feature/) {
			($feature_strand, $feature_start, $feature_end) = getCoordinates ($line);
			push @genes_from_gbk, "$accession\t$feature_strand\t$feature_start\t$feature_end";
		} else {
			#do nothing
		};
	}
	close SEQFILE;

	#calculate upstream for the first gene	
	my @current_gene = split(/\t/, $genes_from_gbk[0]);
	if ($current_gene[1] eq "d") {
		if ($current_gene[2] > 1) {
			my $upstream_end = $current_gene[2] - 1;
			push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t1\t$upstream_end";
		}; 
	} elsif ($current_gene[1] eq "r") {
		my @downstream_gene = split(/\t/, $genes_from_gbk[1]);
		if ($current_gene[3] < $downstream_gene[2] - 1) {
			my $upstream_start = $current_gene[3] +1;
			my $upstream_end = $downstream_gene[2] - 1;
			push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t$upstream_start\t$upstream_end";
#			print "$current_gene[0]\t$current_gene[1]\t$upstream_start\t$upstream_end\n";
		}; 
	}

	#calculate upstreams for the next gene	s
	for (my $i = 1; $i < @genes_from_gbk - 1; $i++) {
		@current_gene = split(/\t/, $genes_from_gbk[$i]);
		my $upstream_start = 0;
		my $upstream_end = 0;
		if ($current_gene[1] eq "d") {
			my @upstream_gene = split(/\t/, $genes_from_gbk[$i-1]);
			$upstream_start = $upstream_gene[3] + 1;
			$upstream_end = $current_gene[2] - 1;
		} elsif ($current_gene[1] eq "r") {
			my @downstream_gene = split(/\t/, $genes_from_gbk[$i+1]);
			$upstream_start = $current_gene[3] + 1;
			$upstream_end = $downstream_gene[2] - 1;
		} else {
			print "ERROR: illegal strand identifier $current_gene[1]";
		}
		if ($upstream_start <= $upstream_end) {
			push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t$upstream_start\t$upstream_end";
		}
	}
	
	#calculate upstream for the last gene	
	@current_gene = split(/\t/, $genes_from_gbk[-1]);
	if ($current_gene[1] eq "d") {
		my @upstream_gene = split(/\t/, $genes_from_gbk[-2]);
		my $upstream_start = $upstream_gene[3] + 1;
		my $upstream_end = $current_gene[2] - 1;
		if ($upstream_start <= $upstream_end) {
			push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t$upstream_start\t$upstream_end";
		}; 
	} elsif ($current_gene[1] eq "r") {
		if ($current_gene[3] < $size) {
			my $upstream_start = $current_gene[3] + 1;
			push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t$upstream_start\t$size";
		}; 
		
	}



	$genes{$accession} = \@genes_from_gbk;
	$upstreams{$accession} = \@upstreams_data;

	$sequences{$accession} = join("", @sequence);
}

sub getCoordinates {
	my ($line) = @_;
	my $feature_strand = "not set";
	my $feature_start = 0;
	my $feature_end = 0;
	$line =~ s/     CDS             //g;
	$line =~ s/\)//g;
	$line =~ s/\>//g;
	$line =~ s/\<//g;
	if ($line =~ /complement/) {
		$feature_strand = "r";
		$line =~ s/complement\(//g;
	} else {
		$feature_strand = "d";
		$line =~ s/\(//g;
	};
	($feature_start, $feature_end) = split(/\.\./, $line);  
	return ($feature_strand, $feature_start, $feature_end);
}

sub process_seq_line {
	my ($line) = @_;
	$line =~ s/[0-9]//g;
	$line =~ s/\s//g;
	return $line;
}

sub get_peak_sequence {
	my ($line) = @_;
	my $ret_val = "";
	my $accession = "";
	my $peak_start = 0;
	my $peak_end = 0;
#	print $line . "\n";
	my @peak_data = split(/\t/, $line);
	$accession = $peak_data[0];
	$accession =~ s/\|.*//g; 
	if (exists $sequences{$accession}) {
		$peak_start = $peak_data[1];
		$peak_end = $peak_data[2];
		if ((!$exclude_peaks_within_genes) || (!is_within_coding_region($accession, $peak_start, $peak_end))){
			if (($peak_end - $peak_start) < $max_fragment_length) {
				$ret_val .= ">" . $accession . "|" . $peak_start . "|" . $peak_end . "\n";
				$ret_val .= substr $sequences{$accession}, $peak_start - 1, $peak_end - $peak_start + 1;
			}
		}
#	} else {
#		print "Accession $accession not found\n";
	}
#	print "\nSeq: ".$ret_val;
	return $ret_val;
}

sub is_within_coding_region {
	my ($accession, $peak_start, $peak_end) = @_;
	my $ret_val = 1;

	#find all upstreams overlapping a peak 
	my @upstreams_data = @{$upstreams{$accession}};
	for (my $i = 0; $i < @upstreams_data; $i++) {
		my @current_upstream = split(/\t/, $upstreams_data[$i]);
		if ((($current_upstream[2]<=$peak_start)&&($current_upstream[3]>=$peak_start))||
			(($current_upstream[2]<=$peak_end)&&($current_upstream[3]>=$peak_end))||
			(($current_upstream[2]>=$peak_start)&&($current_upstream[3]<=$peak_end))) {
			$ret_val = 0;
		}
	}
	return $ret_val;
}

