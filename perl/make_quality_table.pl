use LWP::Simple;
use strict;

my $infile = "quality_mapping.txt";
my $outfile = "quality_mapping_table.tsv";

my $reads_processed = "";
my $reads_aligned = "";
my $reads_failed = "";
my $reads_suppressed = "";
my $data_file = "";

open (OUTFILE, ">$outfile");
print OUTFILE "Data file\tReads processed\tReads aligned\tReads failed to align\tReads with suppressed alignments\n";

open (INFILE, "$infile");

while (my $line = <INFILE>) {
	chomp $line;
	if ($line eq "") { #do nothing
	} elsif ($line =~ /^# reads processed: /){
		$line =~ s/^# reads processed: //;
		$reads_processed = $line;
	} elsif ($line =~ /^# reads with at least one reported alignment: /){
		$line =~ s/^# reads with at least one reported alignment: //;
		$reads_aligned = $line;
	} elsif ($line =~ /^# reads that failed to align: /){
		$line =~ s/^# reads that failed to align: //;
		$reads_failed = $line;
	} elsif ($line =~ /^# reads with alignments suppressed due to -m: /){
		$line =~ s/^# reads with alignments suppressed due to -m: //;
		$reads_suppressed = $line;
	} elsif ($line =~ /^Reported/){
		print OUTFILE "$data_file\t$reads_processed\t$reads_aligned\t$reads_failed\t$reads_suppressed\n";
	} else {
		my @ids = split / /, $line;
		$data_file = $ids[-1];
		$data_file =~ s/.sam$//;
	}
}
close INFILE;

close OUTFILE;

exit(0);
