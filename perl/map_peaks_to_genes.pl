#!/usr/bin/perl

use 5.010;
use strict;
use warnings;

my $infile = "test/NarLp_vs_NarLin_peaks.xls";
my $genome_file = "test/genome_list.txt";
my $outfile = "test/NarLp_vs_NarLin_peaks_genes.txt";

#read least of sequences 
my @genome_seq_list = ();
open (LISTFILE, "$genome_file");
while (my $line = <LISTFILE>) {
	chomp $line;
	push @genome_seq_list,$line;
}
close LISTFILE;

#read list of genes
my @genes_data = ();
my @upstreams_data = ();

foreach my $seqfile (@genome_seq_list){
	if (-e $seqfile){
		print $seqfile."\n";
		processGbkFile($seqfile) ;
	}}

#read infile
open (INFILE, "$infile");
open (OUTFILE, ">$outfile");while (my $line = <INFILE>) {
	chomp $line;
	$line = mapPeaks ($line);
	print OUTFILE $line."\n";
}
close OUTFILE;
close INFILE;
#######################
##### SUBROUTINES #####
#######################
sub processGbkFile {
	my ($seqfile) = @_;
	my $gene_data = "";
	my $gene_flag = 0;
	my $gene_feature = "^     CDS";
	my $gene_qualifier = "                     \/gene=\"";
	my $locus_tag_qualifier = "^                     \/locus_tag=\"";
	my $old_locus_tag_qualifier = "^                     \/old_locus_tag=\"";
	my $feature_start = 0;
	my $feature_end = 0;
	my $feature_strand = "not set";
	my $gene_name = "";
	my $locus_tag = "";
	my $old_locus_tag = " ";
	my $accession = "";
	my $size = 0;
	my @genes_from_gbk = ();

	open (SEQFILE, "$seqfile");
	while (my $line = <SEQFILE>) {
		chomp $line;
		if ($line =~ /^ORIGIN.*/) {
			if ($gene_flag){
				push @genes_from_gbk, "$accession\t$feature_strand\t$feature_start\t$feature_end\t$gene_name\t$locus_tag\t$old_locus_tag";
				$feature_start = 0;
				$feature_end = 0;
				$feature_strand = "not set";
				$gene_name = "";
				$locus_tag = "";
				$old_locus_tag = " ";
				$gene_flag = 0;
			}
			print "Finish parsing $accession sequence\n";
			last;
		} elsif ($line =~ /^LOCUS .*/) {
			my @locus = split(/\s+/, $line);
			$size = $locus[2];
		} elsif ($line =~ /^ACCESSION   /) {
			$line =~ s/^ACCESSION   //g;
			$accession = $line;
		} elsif ($line =~ /$gene_feature/) {
			if ($gene_flag){
				push @genes_from_gbk, "$accession\t$feature_strand\t$feature_start\t$feature_end\t$gene_name\t$locus_tag\t$old_locus_tag";
				$gene_name = "";
				$locus_tag = "";
				$old_locus_tag = " ";
			} else {
				$gene_flag = 1;
			}
			($feature_strand, $feature_start, $feature_end) = getCoordinates ($line);
		} elsif ($line =~ /^\s\s\s\s\s[A-Za-z]/) {
			if ($gene_flag){
				push @genes_from_gbk, "$accession\t$feature_strand\t$feature_start\t$feature_end\t$gene_name\t$locus_tag\t$old_locus_tag";
				$feature_start = 0;
				$feature_end = 0;
				$feature_strand = "not set";
				$gene_name = "";
				$locus_tag = "";
				$old_locus_tag = " ";
				$gene_flag = 0;
			}
		} elsif ($line =~ /$locus_tag_qualifier/) {
			if ($gene_flag){
				$locus_tag = cleanID($line);
			} else {
				#print "Locus tag qualifier found out of CDS feature $line \n";
			};
		} elsif ($line =~ /$old_locus_tag_qualifier/) {
			if ($gene_flag){
				$old_locus_tag = cleanID($line);
			} else {
				#print "Old locus tag qualifier found out of CDS feature $line \n";
			};
		} elsif ($line =~ /$gene_qualifier/) {
			if ($gene_flag){
				$gene_name = cleanID($line);
			} else {
				#print "Gene qualifier found out of CDS feature $line \n";
			};
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
			push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t1\t$upstream_end\t$current_gene[4]\t$current_gene[5]\t$current_gene[6]";
		}; 
	} elsif ($current_gene[1] eq "r") {
		my @downstream_gene = split(/\t/, $genes_from_gbk[1]);
		if ($current_gene[3] < $downstream_gene[2] - 1) {
			my $upstream_start = $current_gene[3] +1;
			my $upstream_end = $downstream_gene[2] - 1;
			push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t$upstream_start\t$upstream_end\t$current_gene[4]\t$current_gene[5]\t$current_gene[6]";
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
			print "ERROR: illegal strand identifier $current_gene[1] for gene $current_gene[5]";
		}
		if ($upstream_start <= $upstream_end) {
			push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t$upstream_start\t$upstream_end\t$current_gene[4]\t$current_gene[5]\t$current_gene[6]";
		}
	}	
	#calculate upstream for the last gene	
	@current_gene = split(/\t/, $genes_from_gbk[-1]);
	if ($current_gene[1] eq "d") {
		my @upstream_gene = split(/\t/, $genes_from_gbk[-2]);
		my $upstream_start = $upstream_gene[3] + 1;
		my $upstream_end = $current_gene[2] - 1;
		if ($upstream_start <= $upstream_end) {
			push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t$upstream_start\t$upstream_end\t$current_gene[4]\t$current_gene[5]\t$current_gene[6]";
		}; 	} elsif ($current_gene[1] eq "r") {
		if ($current_gene[3] < $size) {
			my $upstream_start = $current_gene[3] + 1;
			push @upstreams_data, "$current_gene[0]\t$current_gene[1]\t$upstream_start\t$size\t$current_gene[4]\t$current_gene[5]\t$current_gene[6]";
		}; 
		
	}
	

	
	push @genes_data, @genes_from_gbk;
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

sub cleanID {
	my ($line) = @_;
	$line =~ s/^(.*?)"//;
	$line =~ s/"//;
	return $line;
}

sub mapPeaks {
	my ($line) = @_;
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
	$accession = $peak_data[0];
	$accession =~ s/\.1\|//g; 
	$peak_start = $peak_data[1];
	$peak_end = $peak_data[2];
	$peak_summit = $peak_data[4];
	print "accession $accession summit $peak_summit \n";

#find upstream at peak summit
	my $upstream_flag = 0;
	my $upstream_list="";
	for (my $i = 0; $i < @upstreams_data; $i++) {
		my @current_upstream = split(/\t/, $upstreams_data[$i]);
		if (($current_upstream[0] eq $accession)&&($current_upstream[2]<=$peak_summit)&&($current_upstream[3]>=$peak_summit)) {
			$upstream_list .= "$current_upstream[5] ";
			if ($current_upstream[6] ne " ") {
				$upstream_list.="($current_upstream[6]) ";
			}
			$upstream_flag = 1;
		}
	}
	if ($upstream_flag) {
		$line.="\t$upstream_list";
	} else {
		$line.="\t ";
	}
#find genes at peak summit	
	my $gene_flag = 0;
	my $gene_list="";	for (my $i = 0; $i < @genes_data; $i++) {
		my @current_gene = split(/\t/, $genes_data[$i]);
		if (($current_gene[0] eq $accession)&&($current_gene[2]<=$peak_summit)&&($current_gene[3]>=$peak_summit)) {
			$gene_list.="$current_gene[5] ";
			if ($current_gene[6] ne " ") {
				$gene_list.="($current_gene[6]) ";
			}
			$gene_flag = 1;
		}
	}
	if ($gene_flag) {
		$line.="\t$gene_list";
	} else {
		$line.="\t ";
	}
#find all upstreams overlapping a peak 
	my $upstreams_flag = 0;
	my $upstreams_list="";
	for (my $i = 0; $i < @upstreams_data; $i++) {
		my @current_upstream = split(/\t/, $upstreams_data[$i]);
		if (($current_upstream[0] eq $accession)&&
			((($current_upstream[2]<=$peak_start)&&($current_upstream[3]>=$peak_start))||
			(($current_upstream[2]<=$peak_end)&&($current_upstream[3]>=$peak_end))||
			(($current_upstream[2]>=$peak_start)&&($current_upstream[3]<=$peak_end)))) {
			$upstreams_list .= "$current_upstream[5] ";
			if ($current_upstream[6] ne " ") {
				$upstreams_list.="($current_upstream[6]) ";
			}
			$upstreams_flag = 1;
		}
	}
	if ($upstreams_flag) {
		$line.="\t$upstreams_list";
	} else {
		$line.="\t ";
	}

#find all genes overlapping a peak 
	my $genes_flag = 0;
	my $genes_list="";
	for (my $i = 0; $i < @genes_data; $i++) {
		my @current_gene = split(/\t/, $genes_data[$i]);
		if (($current_gene[0] eq $accession)&&((($current_gene[2]<=$peak_start)&&($current_gene[3]>=$peak_start))||
			 (($current_gene[2]<=$peak_end)&&($current_gene[3]>=$peak_end))||
			 (($current_gene[2]>=$peak_start)&&($current_gene[2]<=$peak_end)))) {
			$genes_list.="$current_gene[5] ";
			if ($current_gene[6] ne " ") {
				$genes_list.="($current_gene[6]) ";
			}
			$genes_flag = 1;
		}
	}
	if ($genes_flag) {
		$line.="\t$genes_list";
	} else {
		$line.="\t ";
	}
	return $line;
}