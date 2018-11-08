#!/public/agis/znjy_group/xianglilingzi/software/miniconda/bin/perl

=head1 NAME

statistics.pl - sequence information statistics

=head1 SYNOPSIS

USAGE: 

  statistics.pl <input_file>

Optional arguments:

    --input_file,-i                     Path to a sequence file
    --reflength,-r                      The length of reference sequence
    --type,-t                           The type of the input file (fastq|fasta, default=auto)
    --head                              Export header
    --help,-h                           Help message

=head1 OPTIONS

B<--input_file,-i>

    A sequence file.

B<--reflength,-r>

    The length of the reference sequence.

B<--type,-t>

    The type of the input file. The default value is auto.

B<--head>

    If you specify this option, the header will be exported.

B<--help,-h>

    This help message

=head1  CONTACT

    Lingzi Xiangli
    xianglilingzi@icloud.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use PerlIO::gzip;
use File::Basename;

my (%options, %length, %GC, %N, @len, $N50, $N90);
my ($sumlength, $sumGC, $sumN, $sum, $coverage)=(0,0,0,0,0);
GetOptions(\%options,
	'input_file|i=s',
	'reflength|r=i',
	'type|t=s',
	'head',
	'help|h'
) || pod2usage();

## display documentation;
if ($options{'help'}){
	pod2usage({-exitval => 0, -verbose => 2, -output => \*STDERR});
}

## make sure everything passed was peachy;
&check_parameters();

## read reference file;
&read_fastq if ($options{type} =~ /fastq/);
&read_fasta if ($options{type} =~ /fasta/);

## statistics;
foreach my $k (keys %length){
	$sumlength += $length{$k};
	push @len, $length{$k};
	$sumGC += $GC{$k};
	$sumN += $N{$k};
}
my $seqnum = @len;
my $GC = $sumGC/($sumlength-$sumN)*100;
my $N = $sumN/$sumlength*100;
if ($options{type} =~ /fasta/){
	my @sortlen = sort {$b<=>$a} @len;
	foreach (@sortlen){
		$sum += $_;
		$N50 = $_ if (!$N50 && $sum > $sumlength*0.5);
		$N90 = $_ if (!$N90 && $sum > $sumlength*0.9);
	}
}
my $coverage = $sumlength/$options{reflength} if ($options{type} =~ /fastq/ && defined $options{reflength});

## output;
my $filename = basename($options{input_file});
if ($options{type} =~ /fasta/){
	printf STDOUT "%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n", 'File name', 'Seq number', 'Total length', 'GC content(%)', 'Gap rate(%)', 'N50 length', 'N90 length' if (defined $options{head});
	printf STDOUT "%15s\t%15d\t%15d\t%15.2f\t%15.2f\t%15d\t%15d\n", $filename, $seqnum, $sumlength, $GC, $N, $N50, $N90;
}
if ($options{type} =~ /fastq/){
	if (defined $options{reflength}){
		printf STDOUT "%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n", 'File name', 'Read number', 'Clean Base', 'GC content(%)', 'Gap rate(%)', 'Coverage' if (defined $options{head});
		printf STDOUT "%15s\t%15d\t%15d\t%15.2f\t%15.2f\t%15.2f\n", $filename, $seqnum, $sumlength, $GC, $N, $coverage;
	}else{
		printf STDOUT "%15s\t%15s\t%15s\t%15s\t%15s\n", 'File name', 'Read number', 'Clean Base', 'GC content(%)', 'Gap rate(%)' if (defined $options{head});
		printf STDOUT "%15s\t%15d\t%15d\t%15.2f\t%15.2f\n", $filename, $seqnum, $sumlength, $GC, $N;
	}
}


#########################
#######Sub Function######
#########################

sub check_parameters {
	die "you must input a --input_file" unless (defined $options{input_file});
	if (!defined $options{type}){
		if ($options{input_file} =~ /\.fq|\.fastq/){
			$options{type}='fastq';
		}
		if ($options{input_file} =~ /\.fa|\.fna|\.fasta/){
			$options{type}='fasta';
		}
	}
	if ($options{input_file} =~ /\.gz$/){
		open IN, "<:gzip", "$options{input_file}" or die $!;
	}else{
		open IN, "$options{input_file}" or die $!;
	}
}

sub read_fasta {
	$/=">"; my $seq = <IN>; $/="\n";
	while (<IN>){
		my $title = $_;
		$/=">";
		(my $seq = <IN>) =~ s/\s+//g;
		chomp $seq;
		$/="\n";
		$length{$title}=length($seq);
		$GC{$title}=( $seq =~ s/[GgCc]/#/g );
		$N{$title}=( $seq =~ s/[Nn]/#/g );
	}
	$/="\n";
	close IN;
}

sub read_fastq {
	while (<IN>){
		my $title = $_;
		(my $seq = <IN>) =~ s/\s+//g;
		$length{$title}=length($seq);
		$GC{$title}=( $seq =~ s/[GgCc]/#/g );
		$N{$title}=( $seq =~ s/[Nn]/#/g );
		<IN>;<IN>;
	}
	close IN;
}
