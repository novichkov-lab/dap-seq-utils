#!/usr/bin/perl

use 5.010;
use strict;
use warnings;
use Excel::Writer::XLSX;
use Cwd 'abs_path';
use File::Basename;
use File::Spec::Functions 'catfile';

# WARNING: This script changes original peak IDs with its owns!

# Paths to files and directories
my $datadir = catfile(dirname(dirname(abs_path($0))),"data");

my $directory_file = "";
my $work_dir = "";
my $motif_dir = "";
my $excel_report_file = "";

if (@ARGV == 4) {
	$directory_file = $ARGV[0];
	$work_dir = $ARGV[1];
	$motif_dir = $ARGV[2];
	$excel_report_file = $ARGV[3];
} else {
	print "Usage: perl generate_report.pl <directory file> <directory with MACS2 output> <directory with MEME output> <report file name>\n";
	print "This script creates final report for DAP-seq pipeline run.\n";
	exit(0);
};

my $genome_library_file = catfile($datadir,"libs.tsv");
my %libs = ();
my $genome_name = "";

# Parameters
my $peaks_filename_suffix = "_peaks.xls";
my $motif_filename_suffix = "_motifs_anr.txt";
my $enrichment_cutoff = 0; # Currently not used
my $pvalue_cutoff = 0; # Currently not used 
my $qvalue_cutoff = 0; # Currently not used
my $disregard_upstreams = 1; #set to 1 if you want all peaks in the XLS file. Set to 0 if you want only peaks with summit located in an upstream region
my $range_left = 500;
my $range_right = 500;
my $orf_overlap = 300;


# Data structures

my %replicons = ();
my %gene_strands = ();
my %operon_first_genes = ();
my %old_locus_tags = ();
my %operon2gene_mappings = ();
my %gene2operon_mappings = ();
my %gene_functions = ();
my %img_ids = ();
my %fitgenomics_ids = ();
my %fitgenomics_functions = ();
my %closest_genes = ();
my %gene_ids = ();

my %samples_directory = ();
my %samples_controls = ();
my %proteins_directory = ();
my %chromosomes = ();
my @genes_data = (); #all genes are here
my @upstreams_data = ();
my %peaks_HoH = (); #data structure that stores all data on peaks in format $peaks_HoH{sample id}{peak id}[nucleotide accession - 0,
# peak start - 1, peak end - 2, peak length - 3, peak summit - 4, peak pile up - 5, -log10(pvalue) - 6, fold enrichment - 7, -log10(qvalue) - 8, peak id - 9, 
# upstream at summit (array ref) - 10, gene at summit (array ref) - 11, all upstreams overlapping the peak (array ref)- 12, all genes overlapping the peak (array ref) - 13, predicted site data - 14]
my %max_peak_number = ();

read_libs($genome_library_file);

if (($work_dir ne "")&&(!(-e $work_dir))) {
	print "Directory $work_dir does not exist!\n";
	exit(1);
} elsif ($work_dir eq "") {
	$work_dir = ".";
};

#open (REPORTFILE, ">$report_file");
#my $report_header = "File\tPeaks\tPeaks in upstreams\tPeaks in genes\tPeaks in plasmids\tPredicted sites\tMotif e-value\n";
#print REPORTFILE $report_header;

#read directory file
open (INFILE, "$directory_file") or die ("File $directory_file not found");
while (my $line = <INFILE>) {
	chomp $line;
	if ($line ne ""){
		my @entry = split(/\t/, $line);
		if ($entry[2]) {
			unless ($genome_name) {
				$genome_name = $entry[1];
			}
			unless ($entry[1] eq $genome_name){
				print "Different genome names in the directory file\n";
				die();
			}
			if (exists $samples_directory{$entry[0]}){
				print "Duplicated entry for sample $entry[0] \n";
				die();
			} else {
				$samples_directory{$entry[0]} = $entry[2];
			}
			if (exists $proteins_directory{$entry[2]}){
				push @{$proteins_directory{$entry[2]}}, $entry[0];
			} else {
				my @arr = ($entry[0]);
				$proteins_directory{$entry[2]} = \@arr;		
			}
		}
		if ($entry[5]){
			$samples_controls{$entry[0]} = $entry[5];
		}
	}
}
close INFILE;

my $gene_function_file = $libs{$genome_name}{"img_data"};
my $operon_file = $libs{$genome_name}{"operons"};
my $fitgenomics_genes_file = $libs{$genome_name}{"fitgenomics_data"};
my $fitgenomics_organism = $libs{$genome_name}{"fitgenomics_name"};
my $replicon_file = $libs{$genome_name}{"replicons"};

my $fitgenomics_url = "http://fit.genomics.lbl.gov/cgi-bin/geneOverview.cgi?orgId=" . $fitgenomics_organism . "&gene=";
my $img_url = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=GeneDetail&page=geneDetail&gene_oid=";

#read list of replicons
load_replicon_data($replicon_file);

#read list of genes from fitgenomics
load_fitgenomics_gene_data($fitgenomics_genes_file);
calculate_upstreams();

#read list of gene functions from IMG
load_img_gene_data($gene_function_file);

#read list of operons
load_operon_data($operon_file);

foreach my $sample (keys %samples_directory) {
	my $file = catfile($work_dir, $sample . "_vs_" . $samples_controls{$sample} . $peaks_filename_suffix);
	if (-e $file){
		print $file . " found\n";
		my $motif_file = catfile($motif_dir, $sample . "_vs_" . $samples_controls{$sample} . $peaks_filename_suffix . $motif_filename_suffix);
		
		#read MACS2 text file (with xls extension), map peaks to genes and calculate statistics
		open (INFILE, "$file") or die ("Unable to open $file");
		my $peaks_count = 0;
		my $peaks_noncoding_count = 0;
		my $peaks_coding_count = 0;
		my $peaks_plasmid = 0;
		while (my $line = <INFILE>) {
			chomp $line;
			if ($line !~ /^\t.*/){
				$line =~ s/\Q$work_dir//g;
				if ($disregard_upstreams) {
					$line = mapPeaks2 ($line, $sample);
				} else {
					$line = mapPeaks ($line, $sample);
				}
			}
		}
		close INFILE;
		get_meme_sites($sample, $motif_file);
	} else {
		print "File $file not found\n";
	}
}

#write excel file
write_excel_report($excel_report_file);


#######################
##### SUBROUTINES #####
#######################


sub read_libs{
	my ($libs_file) = @_;
	open (INFILE, $libs_file) or die ("Unable to open libs file $libs_file");
	while (my $line = <INFILE>) {
		chomp $line;
		my ($lib_name, $lib_path, $size, $fg_file, $fg_name, $img_file, $operons_file, $replicons_file) = split(/\t/, $line);
		$libs{$lib_name}{"path"} = catfile($datadir,$lib_path);
		$libs{$lib_name}{"size"} = $size;
		$libs{$lib_name}{"fitgenomics_data"} = catfile($datadir,$fg_file);
		$libs{$lib_name}{"fitgenomics_name"} = $fg_name;
		$libs{$lib_name}{"img_data"} = catfile($datadir,$img_file);
		$libs{$lib_name}{"operons"} = catfile($datadir,$operons_file);
		$libs{$lib_name}{"replicons"} = catfile($datadir,$replicons_file);
	}
	return 1;
}

sub load_replicon_data {
	my ($replicon_file) = @_;
	open (INFILE, "$replicon_file") or die ("File $replicon_file not found");
	while (my $line = <INFILE>) {
		chomp $line;
		if ($line ne ""){
			my @entry = split(/\t/, $line);
			$replicons{$entry[0]} = $entry[1];
		}
	}
	close INFILE;	
}

sub load_img_gene_data {
	my ($infile) = @_;
	open (INFILE, "$infile") or die ("File $infile not found");
	while (my $line = <INFILE>) {
		chomp $line;
		if ($line ne ""){
			my @entry = split(/\t/, $line);
			$gene_functions{$entry[1]} = $entry[2];
			$img_ids{$entry[1]} = $entry[0];
			$gene_ids{$entry[0]} = $entry[1];
		}
	}
	close INFILE;

}
sub load_fitgenomics_gene_data {
	my ($infile) = @_;
	open (INFILE, "$infile") or die ("File $infile not found");
	my $line = <INFILE>;
	while ($line = <INFILE>) {
		chomp $line;
		if ($line ne ""){
			my @entry = split(/\t/, $line);
			my $strand = "";
			if ($entry[5] eq "+") {
				$strand = "d";
			} elsif ($entry[5] eq "-") {
				$strand = "r";
			} else {
				print "Strand parsing error: $line\n";
			}
			my $gene_dataline = $entry[2] . "\t" . $strand . "\t" . $entry[3] . "\t" . $entry[4] . "\t\t" . $entry[1] . "\t" . $entry[1];
			push @genes_data, $gene_dataline;
			
			$fitgenomics_ids{$entry[1]} = $entry[0];
			$fitgenomics_functions{$entry[1]} = $entry[7];
			$gene_strands{$entry[1]} = $entry[5];
		}
	}
	close INFILE;
}

sub load_operon_data {
	my ($infile) = @_;
	open (INFILE, "$infile") or die ("File $infile not found");
	my $line = <INFILE>;
	while ($line = <INFILE>) {
		chomp $line;
		if ($line ne ""){
			my @entry = split(/\t/, $line);
			my @genes = @entry[5..$#entry];
			if ($entry[3] eq "+") {
				$operon2gene_mappings{$entry[0]} = \@genes;
				foreach my $gene (@genes){
					$gene =~ s/ //g;
					$gene2operon_mappings{$gene} = $entry[0];
				}
			} elsif ($entry[3] eq "-") {
				@genes = reverse @genes;
				$operon_first_genes{$genes[0]} = $entry[0];
				$operon2gene_mappings{$entry[0]} = \@genes;
				foreach my $gene (@genes){
					$gene =~ s/ //g;
					$gene2operon_mappings{$gene} = $entry[0];
				}
			} else {
				print "Error: $entry[3] is not a valid strand designation\n";
			}
		}
	}
	close INFILE;
}


sub calculate_upstreams {

	#calculate upstream for the first gene	
	
	my @current_gene = split(/\t/, $genes_data[0]);
	my $current_accession = $current_gene[0];
	my $upstream_start = 0;
	my $upstream_end = -1;
	my $previous_gene_end = 0;
	
	#calculate upstreams for all genes but the last one
	for (my $i = 0; $i < @genes_data - 1; $i++) {
		@current_gene = split(/\t/, $genes_data[$i]);
		if ($current_gene[0] eq $current_accession) {
			
			if ($current_gene[1] eq "d") {
				$upstream_start = $previous_gene_end + 1;
				$upstream_end = $current_gene[2] - 1;
			} elsif ($current_gene[1] eq "r") {			
				my @downstream_gene = split(/\t/, $genes_data[$i+1]);
				$upstream_start = $current_gene[3] + 1;
				if ($downstream_gene[0] eq $current_accession) {
					$upstream_end = $downstream_gene[2] - 1;
				} else {
					$upstream_end = $replicons{$current_accession};
				
				}
			}
			$previous_gene_end = $current_gene[3];
		} else {
			$previous_gene_end = 0;
			$current_accession = $current_gene[0];
			if ($current_gene[1] eq "d") {
				$upstream_start = $previous_gene_end + 1;
				$upstream_end = $current_gene[2] - 1;
			} elsif ($current_gene[1] eq "r") {
				$upstream_start = $current_gene[3] + 1;
				my @downstream_gene = split(/\t/, $genes_data[$i+1]);
				if ($downstream_gene[0] eq $current_accession) {
					$upstream_end = $downstream_gene[2] - 1;
				} else {
					$upstream_end = $replicons{$current_accession};
				}
			}
			$previous_gene_end = $current_gene[3];
		}
		
		unless (defined $upstream_end){
		}
		if ($upstream_start <= $upstream_end) {
			push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t$upstream_start\t$upstream_end\t$current_gene[4]\t$current_gene[5]\t$current_gene[6]";
		}
	}
	
	#calculate upstream for the last gene	
	@current_gene = split(/\t/, $genes_data[-1]);
	if ($current_gene[0] eq $current_accession) {
		if ($current_gene[1] eq "d") {
			$upstream_start = $previous_gene_end + 1;
			$upstream_end = $current_gene[2] - 1;
		} elsif ($current_gene[1] eq "r") {
			$upstream_start = $current_gene[3] + 1;
			$upstream_end = $replicons{$current_accession};
		}

	} else {
		$current_accession = $current_gene[0];
		if ($current_gene[1] eq "d") {
			$upstream_start = 1;
			$upstream_end = $current_gene[2] - 1;
		} elsif ($current_gene[1] eq "r") {
			$upstream_start = $current_gene[3] + 1;
			$upstream_end = $replicons{$current_accession};
		}
	}
	if ($upstream_start <= $upstream_end) {
		push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t$upstream_start\t$upstream_end\t$current_gene[4]\t$current_gene[5]\t$current_gene[6]";
	}

}

sub cleanID {
	my ($line) = @_;
	$line =~ s/^(.*?)"//;
	$line =~ s/"//;
	return $line;
}

sub get_sample_id {
	my ($peak_id) = @_;
	my @arr = split(/_/, $peak_id);
	pop (@arr);
#	pop (@arr);
	my $ret_val = join("_", @arr);
#	my $ret_val = $arr[0];
	return $ret_val;
}

sub get_peak_number {
	my ($peak_id) = @_;
	my @arr = split(/_/, $peak_id);
	my $ret_val = pop (@arr);
	return $ret_val;
}

# get_gene_function sub uses old locus tag as an argument and returns either gene function (as defined in IMG export data file) or empty string
sub get_gene_function { 
	my ($locus_id) = @_;
	my $ret_val = "Function unknown";
	if (exists $gene_functions{$locus_id}) {
		$ret_val = $gene_functions{$locus_id};
	}
	return $ret_val;
}

# get_gene_strand sub uses old locus tag as an argument and returns either "+" or "-" or empty string
sub get_gene_strand { 
	my ($locus_id) = @_;
	my $ret_val = "";
	if (exists $gene_strands{$locus_id}) {
		$ret_val = $gene_strands{$locus_id};
	}
	return $ret_val;
}

# get_operon_id sub uses old locus tag as an argument and returns either operon id or "no operon"
sub get_operon_id { 
	my ($locus_id) = @_;
	my $ret_val = "no operon";
	if (exists $gene2operon_mappings{$locus_id}) {
		$ret_val = $gene2operon_mappings{$locus_id};
	}
	return $ret_val;
}

# get_fitgenomics_function sub uses old locus tag as an argument and returns either gene function (as defined in fitgenomics genes data file) or empty string
sub get_fitgenomics_function { 
	my ($locus_id) = @_;
	my $ret_val = "Function unknown";
	if (exists $fitgenomics_functions{$locus_id}) {
		$ret_val = $fitgenomics_functions{$locus_id};
	}
	return $ret_val;
}

# get_old_locustag sub uses locus tags as an argument and returns either old locus tag (if defined in GBK files) or empty string
sub get_old_locustag { 
	my ($locus_id) = @_;
	my $ret_val = ""; 
	if (exists $old_locus_tags{$locus_id}) {
		$ret_val = $old_locus_tags{$locus_id};
	}
	$ret_val =~ s/ //g;
#	return $ret_val;
	return $locus_id;
}

# get_img_link sub accepts old locus tags as an argument and returns either link to gene details at IMG (if IMG id defined in %img_ids) or empty string
sub get_img_link { 
	my ($locus_id) = @_;
	my $ret_val = ""; 
	if (exists $img_ids{$locus_id}) {
		$ret_val = $img_url . $img_ids{$locus_id};
	}
	return $ret_val;
}

# get_fitgenomics_link sub accepts old locus tags as an argument and returns either link to gene details at fit.genomics.lbl.gov (if gene id was defined in %fitgenomics_ids) or empty string
sub get_fitgenomics_link { 
	my ($locus_id) = @_;
	my $ret_val = ""; 
	if (exists $fitgenomics_ids{$locus_id}) {
		$ret_val = $fitgenomics_url . $fitgenomics_ids{$locus_id};
	}
	return $ret_val;
}


sub mapPeaks {
	my ($line, $sample_id) = @_;
	my $accession = "";
	my $peak_start = 0;
	my $peak_end = 0;
	my $peak_summit = 0;
	my $header = "chr	start	end	length	abs_summit	pileup	-log10(pvalue)	fold_enrichment	-log10(qvalue)	name";
	if (($line =~ /^#.*/)||($line eq "")) {
		return $line;
	}
	if ($line eq $header ) {
		$line = $line."\tupstream(s) at summit\tgene(s) at summit\tupstream(s) overlapping peak\tgene(s) overlapping peak";
		return $line;
	}
	my @peak_data = split(/\t/, $line);
	if (($peak_data[7] < $enrichment_cutoff) || ($peak_data[6] < $pvalue_cutoff) || ($peak_data[8] < $qvalue_cutoff)){
		return $line;
	}
	my $peak_number = $peak_data[9];
	$peak_number =~ s/.*_//;
	$max_peak_number{$sample_id} = $peak_number;
	$peak_data[9] = $sample_id . "_" . $peak_number;

#	print $peak_data[9] . "\n";
#	my $sample_id = get_sample_id($peak_data[9]);
	$accession = $peak_data[0];
	$accession =~ s/\|.*//g; 
	$peak_data[0] = $accession;
	$peaks_HoH{$sample_id}{$peak_data[9]} = \@peak_data;
	$peak_start = $peak_data[1];
	$peak_end = $peak_data[2];
	$peak_summit = $peak_data[4];
#	print "accession $accession summit $peak_summit \n";

#find upstream at peak summit
	my $upstream_flag = 0;
#	my @upstream_list = ();
	my @upstreams_at_summit = ();
	$line.= "\t";
	for (my $i = 0; $i < @upstreams_data; $i++) {
		my @current_upstream = split(/\t/, $upstreams_data[$i]);
		if (($current_upstream[0] eq $accession)&&($current_upstream[2]<=$peak_summit)&&($current_upstream[3]>=$peak_summit)) {
			my $upstream_line = "$current_upstream[5]\[";
			push @upstreams_at_summit, $current_upstream[5];
			my @arr = ();
			if (exists $gene_functions{$current_upstream[5]}){
				push @arr, $gene_functions{$current_upstream[5]};
			}
			if ($current_upstream[6] ne " ") {
				push @arr, $current_upstream[6];
				if (exists $gene_functions{$current_upstream[6]}){
					push @arr, $gene_functions{$current_upstream[6]};
				}
			}
			if (@arr) {
				$upstream_line .= join(", ", @arr);
				$upstream_line .= "\] ";
			} else {
				$upstream_line .= "Unknown function\] ";
			}
			$upstream_flag = 1;
			$line.= $upstream_line . " ";
#			push @upstream_list, $upstream_line;
		}
	}
	if ($upstream_flag) {
		push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, \@upstreams_at_summit); #join ("\n", @upstream_list));
	} else {
		$line.=" ";
		push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, "");
	}

#find genes at peak summit	
	my $gene_flag = 0;
	my $gene_list="";
	my $closest_gene = "";
	my $smallest_distance = 1000000;
	my @genes_at_summit = ();
	for (my $i = 0; $i < @genes_data; $i++) {
		my @current_gene = split(/\t/, $genes_data[$i]);
		my $include_gene_flag = 0;
		my $gene_strand = $current_gene[1];
		my $gene_left = $current_gene[2];
		my $gene_right = $current_gene[3];
		if ($current_gene[0] eq $accession) {
			if ($gene_strand eq "d") {
				my $distance = $peak_summit - $gene_left;
				if ($distance < 0) {
					$distance = 0 - $distance;
				};
				if ($distance < $smallest_distance) {
					$smallest_distance = $distance;
					$closest_gene = $current_gene[5];
				}

				my $right_limit = $gene_left + $orf_overlap;
				if ($right_limit > $gene_right) {
					$right_limit = $gene_right;
				}
				if (($gene_left<=$peak_summit)&&($right_limit>=$peak_summit)) {
					$include_gene_flag = 1;

				}
			
			} elsif ($gene_strand eq "r") {
				my $distance = $peak_summit - $gene_right;
				if ($distance < 0) {
					$distance = 0 - $distance;
				};
				if ($distance < $smallest_distance) {
					$smallest_distance = $distance;
					$closest_gene = $current_gene[5];
				}

				my $left_limit = $gene_right - $orf_overlap;
				if ($left_limit < $gene_left) {
					$left_limit = $gene_left;
				}
				if (($left_limit<=$peak_summit)&&($gene_right>=$peak_summit)) {
					$include_gene_flag = 1;
				}
			}
		}
		
		if ($include_gene_flag) {
			push @genes_at_summit, $current_gene[5];
			$gene_list.="$current_gene[5]\[";
			my @arr = ();
			if (exists $gene_functions{$current_gene[5]}){
				push @arr, $gene_functions{$current_gene[5]};
			}
			if ($current_gene[6] ne " ") {
				push @arr, $current_gene[6];
				if (exists $gene_functions{$current_gene[6]}){
					push @arr, $gene_functions{$current_gene[6]};
				}
			}
			if (@arr) {
				$gene_list .= join(", ", @arr);
				$gene_list .= "\] ";
			} else {
				$gene_list .= "Unknown function\] ";
			}

			$gene_flag = 1;
		}
	}
	if ($gene_flag) {
		$line.="\t$gene_list";
		push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, \@genes_at_summit);
	} else {
		$line.="\t ";
		push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, "");
	}
	$closest_genes{$sample_id}{$peak_data[9]} = $closest_gene;

#find all upstreams overlapping a peak 
	my $upstreams_flag = 0;
	my $upstreams_list = "";
	my @upstreams_at_peak = ();
	for (my $i = 0; $i < @upstreams_data; $i++) {
		my @current_upstream = split(/\t/, $upstreams_data[$i]);
		if (($current_upstream[0] eq $accession)&&
			((($current_upstream[2]<=$peak_start)&&($current_upstream[3]>=$peak_start))||
			(($current_upstream[2]<=$peak_end)&&($current_upstream[3]>=$peak_end))||
			(($current_upstream[2]>=$peak_start)&&($current_upstream[3]<=$peak_end)))) {
			$upstreams_list .= "$current_upstream[5] ";
			push @upstreams_at_peak, $current_upstream[5];
			if ($current_upstream[6] ne " ") {
				$upstreams_list.="($current_upstream[6]) ";
			}
			$upstreams_flag = 1;
		}
	}
	if ($upstreams_flag) {
		$line.="\t$upstreams_list";
		push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, \@upstreams_at_peak);
	} else {
		$line.="\t ";
		push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, "");
	}


#find all genes overlapping a peak 
	my $genes_flag = 0;
	my $genes_list="";
	my @genes_at_peak = ();
	for (my $i = 0; $i < @genes_data; $i++) {
		my @current_gene = split(/\t/, $genes_data[$i]);
		if (($current_gene[0] eq $accession)&&((($current_gene[2]<=$peak_start)&&($current_gene[3]>=$peak_start))||
			 (($current_gene[2]<=$peak_end)&&($current_gene[3]>=$peak_end))||
			 (($current_gene[2]>=$peak_start)&&($current_gene[2]<=$peak_end)))) {
			$genes_list.="$current_gene[5] ";
			push @genes_at_peak,$current_gene[5];
			if ($current_gene[6] ne " ") {
				$genes_list.="($current_gene[6]) ";
			}
			$genes_flag = 1;
		}
	}
	if ($genes_flag) {
		$line.="\t$genes_list";
		
		push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, \@genes_at_peak);
	} else {
		$line.="\t ";
		push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, "");
	}
	
	push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, "");
	return $line;
}

# mapping procedure for finding genes in a given range around peak summit 
sub mapPeaks2 {
	my ($line, $sample_id) = @_;
	my $accession = "";
	my $peak_summit = 0;
	my $header = "chr	start	end	length	abs_summit	pileup	-log10(pvalue)	fold_enrichment	-log10(qvalue)	name";
	if (($line =~ /^#.*/)||($line eq "")) {
		return $line;
	}
	if ($line eq $header ) {
		$line = $line."\tupstream(s) at summit\tgene(s) at summit\tupstream(s) overlapping peak\tgene(s) overlapping peak";
		return $line;
	}
	my @peak_data = split(/\t/, $line);
	if (($peak_data[7] < $enrichment_cutoff) || ($peak_data[6] < $pvalue_cutoff) || ($peak_data[8] < $qvalue_cutoff)){
		return $line;
	}
	my $peak_number = $peak_data[9];
	$peak_number =~ s/.*_//;
	$max_peak_number{$sample_id} = $peak_number;
	$peak_data[9] = $sample_id . "_" . $peak_number; 
#	my $sample_id = get_sample_id($peak_data[9]);	
	$accession = $peak_data[0];
	#$accession =~ s/\.1\|//g; 
	$accession =~ s/\|.*//g; 
	$peak_data[0] = $accession;
	$peaks_HoH{$sample_id}{$peak_data[9]} = \@peak_data;
	$peak_summit = $peak_data[4];
	my $left_limit = $peak_summit - $range_left;
	my $right_limit = $peak_summit + $range_right;

	push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, "");

#find genes at peak summit	
	my $gene_flag = 0;
	my $closest_gene = "";
	my $smallest_distance = 1000000;
	my @genes_at_summit = ();
	for (my $i = 0; $i < @genes_data; $i++) {
#		print $genes_data[$i] . "\n" ;
		my @current_gene = split(/\t/, $genes_data[$i]);
		my $include_gene_flag = 0;
		my $gene_strand = $current_gene[1];
		my $gene_left = $current_gene[2];
		my $gene_right = $current_gene[3];
		if ($current_gene[0] eq $accession) {
			if ($gene_strand eq "d") {
				if (($left_limit<=$gene_left)&&($right_limit>=$gene_left) && ($peak_summit < ($gene_left + $orf_overlap - 1))) {
					$include_gene_flag = 1;
					my $distance = $peak_summit - $gene_left;
					if ($distance < 0) {
						$distance = 0 - $distance;
					};
					if ($distance < $smallest_distance) {
						$smallest_distance = $distance;
						$closest_gene = $current_gene[5];
					}
					
				}
			
			} elsif ($gene_strand eq "r") {
				if (($left_limit<=$gene_right)&&($right_limit>=$gene_right) && ($peak_summit > ($gene_left - $orf_overlap + 1))) {
					$include_gene_flag = 1;
					my $distance = $peak_summit - $gene_right;
					if ($distance < 0) {
						$distance = 0 - $distance;
					};
					if ($distance < $smallest_distance) {
						$smallest_distance = $distance;
						$closest_gene = $current_gene[5];
					}

				}
			}
		}
		
		if ($include_gene_flag) {
			push @genes_at_summit, $current_gene[5];
			$gene_flag = 1;
		}
		
		#break the cycle
#		if (($gene_left > $peak_summit + $range_right) && ($gene_right > $peak_summit + $range_right)) {
#			last;
#		}
	}
	if ($gene_flag) {
#		print join (",", @genes_at_summit);
		push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, \@genes_at_summit);
	} else {
		push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, "");
	}
	$closest_genes{$sample_id}{$peak_data[9]} = $closest_gene;
#skip finding all upstreams overlapping a peak 
	push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, "");


#skip finding all genes overlapping a peak 
	push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, "");
	
	push (@{$peaks_HoH{$sample_id}{$peak_data[9]}}, "");

	return $line;
}

sub get_meme_sites{
	my ($sample, $motif_file) = @_;
	my $site_section_start = "Motif 1 sites sorted by position p-value";
	my $site_section_end = "Motif 1 block diagrams";
	my $flag = 0;
	
#	if (-e $motif_file){
#		print $motif_file." found\n";
#	}
	open (INFILE, "$motif_file") or return;
	while (my $line = <INFILE>) {
		chomp $line;
		if ($flag) {
			if ($line =~ /$site_section_end/){
				$flag = 0;
			} else {
				my @fields = split(/\s+/, $line);
				if ((scalar @fields) == 7){
					my ($accession, $start, $end) = split(/\|/, $fields[0]);
					#$accession.= ".1";
					if (exists $peaks_HoH{$sample}) {
						while (my ($peak_id, $peak_data) = each %{ $peaks_HoH{$sample} } ) { 
							my $peak_accession = ${$peak_data}[0];
							#$peak_accession =~ s/\.1\|//;
							if (($peak_accession eq $accession) && (${$peak_data}[1] eq $start)) {
								my $site_start = $start + $fields[2];
								${$peaks_HoH{$sample}{$peak_id}}[14] = $site_start . ":". $fields[5];
#								print $site_start . ":". $fields[5] . "\n";
							} else {
#								print ${$peak_data}[0] . ":" . $accession . ",". ${$peak_data}[1] . ":" . $start  . "\n";
							}
						}
					} else {
						print "Sample ID $sample unknown\n";
					}
				}
			}
		} elsif ($line =~ /$site_section_start/){
			$flag = 1;
		}
	}
	print "Processing of $motif_file completed\n";
	close INFILE;
	
}


sub write_excel_comparison_table{
	my ($workbook, $name, $samples_list_ref) = @_;
	my @samples_list = @{$samples_list_ref};
	
	my $worksheet = $workbook->add_worksheet( $name );
	
	$worksheet->set_column(0, 0, 30);

	if (@samples_list) {
		$worksheet->write(0, 0, "Accession:Start..End");
	} else {
		$worksheet->write(0, 0, "No peaks found in all samples");
		return;
	}
	my $i = 1; #here $i is number of a column
	# define positions of samples on the sheet
	my %sample_ids = ();
	foreach my $sample_id (@samples_list) {
		$sample_ids{$sample_id} = $i;
		$sample_id =~ s/.*\///;
		$worksheet->write(0, $i, $sample_id);
		$i++;
	}
	# arrange peak starts in ascending order
	my %peak_starts = ();
	foreach my $sample_id (@samples_list){
#		print $sample_id . "\n";
		if (scalar (keys %{$peaks_HoH{$sample_id}}) == 0){
			$worksheet->write(1, $sample_ids{$sample_id}, "No peaks found");
			next;
		}
		foreach my $peak_id (keys %{$peaks_HoH{$sample_id}}){
#			print $peak_id . "\n";
			my $peak_start = ${$peaks_HoH{$sample_id}{$peak_id}}[1];
			my $accession = ${$peaks_HoH{$sample_id}{$peak_id}}[0];
			if (exists $peak_starts{$accession}{$peak_start}){
				push @{$peak_starts{$accession}{$peak_start}}, $peak_id;
			} else {
				my @arr = ($peak_id);
				$peak_starts{$accession}{$peak_start} = \@arr;
			}
#			print $peak_start . "\n";
		}
	}
	if (! keys %peak_starts) {
		$worksheet->write(0, 0, "No peaks found");
		return;
	}
	my %peaks_accounted = ();
	my %first_column = ();
	$i = 1; #here $i is a number of row
	my @peaks_cluster = ();
	my $last_accession = "";
	my $cluster_start = -1;
	my $cluster_end = -1;
	foreach my $accession (sort keys %peak_starts){
		$last_accession = $accession;
		@peaks_cluster = ();
		$cluster_start = -1;
		$cluster_end = -1;
#		my @starts_list_global = sort keys %{$peak_starts{$accession}};
		foreach my $position (sort {$a <=> $b} keys %{$peak_starts{$accession}}){
#			print $position."\n";
			foreach my $peak_id (@{$peak_starts{$accession}{$position}}){
#				print "iterate over peaks:" . $peak_id . "\n";
				unless (exists $peaks_accounted{$peak_id}){
					my $sample_id = get_sample_id($peak_id);
					my $peak_start = ${$peaks_HoH{$sample_id}{$peak_id}}[1];
					my $peak_end = ${$peaks_HoH{$sample_id}{$peak_id}}[2];
					if ($cluster_end < 0) {
						# very first cluster
						$cluster_start = $peak_start;
						$cluster_end = $peak_end;
						push @peaks_cluster, $peak_id;
					} elsif ($peak_start > $cluster_end) {
						# write current cluster, 
						my %cells=(); 
						foreach my $cluster_member (@peaks_cluster){
							my $cluster_member_sample_id = get_sample_id($cluster_member);
							if (exists $cells{$cluster_member_sample_id}) {
								$cells{$cluster_member_sample_id} = $cells{$cluster_member_sample_id} + 1;
							} else {
								$cells{$cluster_member_sample_id} = $i;
							}
							$worksheet->write ($cells{$cluster_member_sample_id}, $sample_ids{$cluster_member_sample_id}, get_peak_number($cluster_member));
							$first_column{$cells{$cluster_member_sample_id}} = $accession . ":" . $cluster_start . ".." . $cluster_end; #${$peaks_HoH{$cluster_member_sample_id}{$cluster_member}}[1] ;
						}
						foreach my $sample_id2 (keys %cells){
							if ($cells{$sample_id2} > $i) {
								$i = $cells{$sample_id2};
							}
						}
						$i++;
						foreach my $sample_id2 (keys %sample_ids) {
							$worksheet->write ($i, $sample_ids{$sample_id2}, "***");
						}
						$i++;
						
						# start new cluster
						@peaks_cluster = ($peak_id);
						$cluster_start = $peak_start;
						$cluster_end = $peak_end;
					} elsif (($peak_start == $cluster_start)||(($peak_start > $cluster_start) && ($peak_start <= $cluster_end))){
						# add the peak to current cluster
						push @peaks_cluster, $peak_id;
						$cluster_end = $peak_end;
					}
					$peaks_accounted{$peak_id} = 1;
				}
				
			}
		}
		# write the last cluster for this accession
		if (@peaks_cluster) {
			my %cells=(); 
			foreach my $cluster_member (@peaks_cluster){
				my $cluster_member_sample_id = get_sample_id($cluster_member);
				if (exists $cells{$cluster_member_sample_id}) {
					$cells{$cluster_member_sample_id} = $cells{$cluster_member_sample_id} + 1;
				} else {
					$cells{$cluster_member_sample_id} = $i;
				}
				$worksheet->write ($cells{$cluster_member_sample_id}, $sample_ids{$cluster_member_sample_id}, get_peak_number($cluster_member));
				$first_column{$cells{$cluster_member_sample_id}} = $accession . ":" . $cluster_start . ".." . $cluster_end;#. ${$peaks_HoH{$cluster_member_sample_id}{$cluster_member}}[1] ;
			}
			foreach my $sample_id2 (keys %cells){
				if ($cells{$sample_id2} > $i) {
					$i = $cells{$sample_id2};
				}
			}
			$i++;
			foreach my $sample_id2 (keys %sample_ids) {
				$worksheet->write ($i, $sample_ids{$sample_id2}, "***");
			}
			$i++;
			@peaks_cluster = ();
			$cluster_start = -1;
			$cluster_end = -1;
		}
	}
	

	# write the first column
	foreach my $row_number (keys %first_column){
		$worksheet->write ($row_number, 0, $first_column{$row_number});
	}
	

}

sub write_excel_allpeaks_long_table {
	
	my ($workbook, $samples_list_ref) = @_;

	my $worksheet = $workbook->add_worksheet( 'All_peaks');

	my @sheet1_header = (
			"Peak", # 0
			"Accession", # 1
			"Start", # 2
			"End",  # 3
			"Length", # 4
			"Summit", # 5
			"Pile up", # 6
			"\'-log10(pvalue)", # 7
			"Fold enrichment", # 8
			"\'-log10(qvalue)", # 9
			"MACSD2 peak id", # 10
			"Upstream at summit: locus tag", # 11
			"Upstream at summit: old locus tag", # 12
			"Upstream at summit: strand", # 13
			"Upstream at summit: operon", # 14
			"Upstream at summit: function", # 15
			"Upstream at summit: IMG link", # 16
			"Upstream at summit: fitgenomics function", # 17
			"Upstream at summit: fitgenomics link", # 18
			"Gene at summit: locus tag", # 19
			"Gene at summit: old locus tag", # 20
			"Gene at summit: strand", # 21
			"Gene at summit: strand", # 22
			"Gene at summit: function", # 23
			"Gene at summit: IMG link", # 24
			"Gene at summit: fitgenomics function", # 25
			"Gene at summit: fitgenomics link", # 26
			"All upstreams: locus tag", # 27
			"All upstreams: old locus tag", # 28
			"All upstreams: strand", # 29
			"All upstreams: operon", # 30
			"All upstreams: function", # 31
			"All upstreams: IMG link", # 32
			"All upstreams: fitgenomics function", # 33
			"All upstreams: fitgenomics link", # 34
			"All genes: locus tag", # 35
			"All genes: old locus tag", # 36
			"All genes: strand", # 37
			"All genes: operon", # 38
			"All genes: function", # 39
			"All genes: IMG link", # 40
			"All genes: fitgenomics function", # 41
			"All genes: fitgenomics link", # 42
			"Predicted site"# 43
			);
	my $sheet1_header_ref = \@sheet1_header;
	$worksheet->write_row(0, 0, $sheet1_header_ref);

	$worksheet->set_column(0, 3, 10);
	my $format = $workbook->add_format();
	$format->set_bold();

	my $i = 1; #here $i is number of a row

	foreach my $sample_id (@{$samples_list_ref}){
		$worksheet->write($i, 0, $sample_id,  $format);
		$i++;
#		foreach my $peak_id (keys %{$peaks_HoH{$sample_id}}){
		for (my $peak_number = 1; $peak_number <= $max_peak_number{$sample_id}; $peak_number++) {
#			my $peak_id = $sample_id . "_peak_" . $peak_number;
			my $peak_id = $sample_id . "_" . $peak_number;
			unless (exists $peaks_HoH{$sample_id}{$peak_id}) {
				next;
			}
			my $row_span = 1;
			my $k = $i;
			my $j = 10;  
			if (${$peaks_HoH{$sample_id}{$peak_id}}[$j]) {
				if ($row_span < scalar(@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]})) {
					$row_span = scalar(@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]});
				}
				foreach my $locus_tag (@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]}) {
					$worksheet->write($k, 11, $locus_tag);
					$worksheet->write($k, 12, get_old_locustag($locus_tag));
					$worksheet->write($k, 13, get_gene_strand(get_old_locustag($locus_tag)));
					$worksheet->write($k, 14, get_operon_id(get_old_locustag($locus_tag)));
					$worksheet->write($k, 15, get_gene_function(get_old_locustag($locus_tag)));
					$worksheet->write($k, 16, get_img_link(get_old_locustag($locus_tag)));
					$worksheet->write($k, 17, get_fitgenomics_function(get_old_locustag($locus_tag)));
					$worksheet->write($k, 18, get_fitgenomics_link(get_old_locustag($locus_tag)));
					$k++;
				}
			}
			$j++;
			$k = $i;
			if (${$peaks_HoH{$sample_id}{$peak_id}}[$j]) {
				if ($row_span < scalar(@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]})) {
					$row_span = scalar(@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]});
				}
				foreach my $locus_tag (@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]}) {
					$worksheet->write($k, 19, $locus_tag);
					$worksheet->write($k, 20, get_old_locustag($locus_tag));
					$worksheet->write($k, 21, get_gene_strand(get_old_locustag($locus_tag)));
					$worksheet->write($k, 22, get_operon_id(get_old_locustag($locus_tag)));
					$worksheet->write($k, 23, get_gene_function(get_old_locustag($locus_tag)));
					$worksheet->write($k, 24, get_img_link(get_old_locustag($locus_tag)));
					$worksheet->write($k, 25, get_fitgenomics_function(get_old_locustag($locus_tag)));
					$worksheet->write($k, 26, get_fitgenomics_link(get_old_locustag($locus_tag)));
					$k++;
				}
			}
			$j++;
			$k = $i;
			if (${$peaks_HoH{$sample_id}{$peak_id}}[$j]) {
				if ($row_span < scalar(@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]})) {
					$row_span = scalar(@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]});
				}
				foreach my $locus_tag (@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]}) {
					$worksheet->write($k, 27, $locus_tag);
					$worksheet->write($k, 28, get_old_locustag($locus_tag));
					$worksheet->write($k, 29, get_gene_strand(get_old_locustag($locus_tag)));
					$worksheet->write($k, 30, get_operon_id(get_old_locustag($locus_tag)));
					$worksheet->write($k, 31, get_gene_function(get_old_locustag($locus_tag)));
					$worksheet->write($k, 32, get_img_link(get_old_locustag($locus_tag)));
					$worksheet->write($k, 33, get_fitgenomics_function(get_old_locustag($locus_tag)));
					$worksheet->write($k, 34, get_fitgenomics_link(get_old_locustag($locus_tag)));
					$k++;
				}
			}
			$j++;
			$k = $i;
			if (${$peaks_HoH{$sample_id}{$peak_id}}[$j]) {
				if ($row_span < scalar(@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]})) {
					$row_span = scalar(@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]});
				}
				foreach my $locus_tag (@{${$peaks_HoH{$sample_id}{$peak_id}}[$j]}) {
					$worksheet->write($k, 35, $locus_tag);
					$worksheet->write($k, 36, get_old_locustag($locus_tag));
					$worksheet->write($k, 37, get_gene_strand(get_old_locustag($locus_tag)));
					$worksheet->write($k, 38, get_operon_id(get_old_locustag($locus_tag)));
					$worksheet->write($k, 39, get_gene_function(get_old_locustag($locus_tag)));
					$worksheet->write($k, 40, get_img_link(get_old_locustag($locus_tag)));
					$worksheet->write($k, 41, get_fitgenomics_function(get_old_locustag($locus_tag)));
					$worksheet->write($k, 42, get_fitgenomics_link(get_old_locustag($locus_tag)));
					$k++;
				}
			}
			$j++;
		
			$row_span += $i;
			while ($i < $row_span){
				$j = 0;
				$worksheet->write($i, $j, $peak_id);
#				$worksheet->write($i, $j, "peak_" . $peak_number);
				while ($j < 10){
					$j++;
					$worksheet->write($i, $j, ${$peaks_HoH{$sample_id}{$peak_id}}[$j-1]);
				}
				$worksheet->write($i, 43, ${$peaks_HoH{$sample_id}{$peak_id}}[14]);
				$i++;
			}

		}
	}


}

sub write_excel_allpeaks_table {
	
	my ($workbook, $samples_list_ref) = @_;

	my $worksheet = $workbook->add_worksheet('All_peaks');

	my @sheet1_header = ( #						Column index
			"Peak", # 							0
			"Accession", # 						1
			"Start", # 							2
			"End",  # 							3
			"Length", # 						4
			"Summit", # 						5
			"Pile up", # 						6
			"\'-log10(pvalue)", # 				7
			"Fold enrichment", # 				8
			"\'-log10(qvalue)", # 				9
			"MACSD2 peak id", # 				10
			"Peak upstream/downstream of:", # 	11
			"Strand", # 						12
			"operon number", # 					13
			"peak closest to:", #				14
			"Function", #						15
			"Predicted site" #					16
			);
	my $sheet1_header_ref = \@sheet1_header;
	$worksheet->write_row(0, 0, $sheet1_header_ref);

	$worksheet->set_column(0, 3, 10);
	my $format = $workbook->add_format();
	$format->set_bold();

	my $i = 1; #here $i is number of a row

	foreach my $sample_id (@{$samples_list_ref}){
#		print $sample_id . "\n";
		if (exists $gene_ids{$samples_directory{$sample_id}}){
			$worksheet->write($i, 0, $sample_id . " (" . $samples_directory{$sample_id} . "-". $gene_ids{$samples_directory{$sample_id}} .")",  $format);
		} else {
			$worksheet->write($i, 0, $sample_id . " (" . $samples_directory{$sample_id} . ")",  $format);
		}
		$i++;
#		foreach my $peak_id (sort {$a<=>$b} keys %{$peaks_HoH{$sample_id}}){
		for (my $peak_number = 1; $peak_number <= $max_peak_number{$sample_id}; $peak_number++) {
#			my $peak_id = $sample_id . "_peak_" . $peak_number;
			my $peak_id = $sample_id . "_" . $peak_number;
#			print $peak_id . "\n";
			unless (defined $peaks_HoH{$sample_id}{$peak_id}) {
				next;
			}
			my $row_span = 0;
			my %peak_genes_list = ();
			my %closest_genes_list = ();

			if (${$peaks_HoH{$sample_id}{$peak_id}}[10]) {
				foreach my $locus_tag (@{${$peaks_HoH{$sample_id}{$peak_id}}[10]}) {
					my $old_lt = get_old_locustag($locus_tag);
#					print "\"" . $locus_tag . "\":\t" . $old_lt . "\"\n";
					if (exists $operon_first_genes{$old_lt}) {
						foreach my $operon_member (@{$operon2gene_mappings{$operon_first_genes{$old_lt}}}){
							$peak_genes_list{$operon_member} = 1;
						}
					} elsif (exists $gene2operon_mappings{$old_lt}) {
						my $flag = 0;
						foreach my $operon_member (@{$operon2gene_mappings{$gene2operon_mappings{$old_lt}}}) {
							if ($operon_member eq $old_lt) {
								$flag = 1;
								$peak_genes_list{$operon_member} = 1;
							} elsif ($flag) {
								$peak_genes_list{$operon_member} = 1;
							}
						}
					} else {
						$peak_genes_list{$old_lt} = 1;
					}
				}
			}
			if (${$peaks_HoH{$sample_id}{$peak_id}}[11]) {
				foreach my $locus_tag (@{${$peaks_HoH{$sample_id}{$peak_id}}[11]}) {
					#my $closest_flag = 0;
					#if ($closest_genes{$sample_id}{$peak_id} eq $locus_tag) {
						#$closest_flag = 1;
					#}
					my $old_lt = get_old_locustag($locus_tag);
#					print "\"" . $locus_tag . "\":\t" . $old_lt . "\"\n";
					if (exists $gene2operon_mappings{$old_lt}) {
						my $flag = 0;
						foreach my $operon_member (@{$operon2gene_mappings{$gene2operon_mappings{$old_lt}}}) {
							if ($operon_member eq $old_lt) {
								$flag = 1;
								$peak_genes_list{$operon_member} = 1;
								#if ($closest_flag) {
									#$closest_genes_list{$operon_member} = $operon_member;
								#} else {
									#$closest_genes_list{$operon_member} = "";
								#}
							} elsif ($flag) {
								$peak_genes_list{$operon_member} = 1;
								#if ($closest_flag) {
									#$closest_genes_list{$operon_member} = $operon_member;
								#} else {
									#$closest_genes_list{$operon_member} = "";
								#}
							}
						}
					} else {
						$peak_genes_list{$old_lt} = 1;
						#if ($closest_flag) {
							#$closest_genes_list{$old_lt} = $old_lt;
						#} else {
							#$closest_genes_list{$old_lt} = "";
						#}
					}
				}
			}
			foreach my $old_locus_tag (keys %peak_genes_list) {
				if ($closest_genes{$sample_id}{$peak_id} eq $old_locus_tag) {
					$closest_genes_list{$old_locus_tag} = $old_locus_tag;
				} else {
					$closest_genes_list{$old_locus_tag} = "";
				}
			}
			my $k = $i;
			foreach my $old_locus_tag (sort keys %peak_genes_list) {
				if ($old_locus_tag ne "") {
					#$worksheet->write($k, 11, $old_locus_tag);
					$worksheet->write($k, 11, $img_ids{$old_locus_tag});
					$worksheet->write($k, 12, get_gene_strand($old_locus_tag));
					$worksheet->write($k, 13, get_operon_id($old_locus_tag));
					#$worksheet->write($k, 14, $closest_genes_list{$old_locus_tag});
					$worksheet->write($k, 14, $img_ids{$closest_genes_list{$old_locus_tag}});
					if ($closest_genes_list{$old_locus_tag} ne ""){
						$worksheet->write($k, 15, get_gene_function($closest_genes_list{$old_locus_tag}));
					} else {
						$worksheet->write($k, 15, "");
					}
					$k++;
					$row_span++;
				}
			}
			if (!$row_span) {
				$row_span++;
			}

			$row_span += $i;
			while ($i < $row_span){
				my $j = 0;
				$worksheet->write($i, $j, $peak_id);
				while ($j < 10){
					$j++;
					$worksheet->write($i, $j, ${$peaks_HoH{$sample_id}{$peak_id}}[$j-1]);
				}
				$worksheet->write($i, 16, ${$peaks_HoH{$sample_id}{$peak_id}}[14]);
				$i++;
			}

		}
	}


}

sub write_excel_report{
	my ($outfile) = @_;
	print "Outfile: ".$outfile."\n";
	my $workbook = Excel::Writer::XLSX->new("$outfile");
	
	my @samples_list = sort keys %peaks_HoH;

	write_excel_allpeaks_table($workbook, \@samples_list);
#	write_excel_allpeaks_long_table($workbook, \@samples_list);

#	write_excel_comparison_table($workbook, 'Peak_comparison', \@samples_list);

	#Create protein-specific pages
	print "Start creating protein pages\n";
	foreach my $protein (sort keys %proteins_directory){
#		print $protein . "\n";
		
		@samples_list = ();
		foreach my $sample (@{$proteins_directory{$protein}}){
			foreach my $sample_id (keys %peaks_HoH) {
				if ($sample_id eq $sample) {#if ($sample_id =~ /\Q$sample/) {
					push @samples_list, $sample_id;
#					last;
				}
			}
		}
#		print "Send to sub:" . join (",", @samples_list) . "\n";
		write_excel_comparison_table($workbook, $protein, \@samples_list);
	}
	
	
	$workbook->close;
}
