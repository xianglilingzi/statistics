#!/public/agis/znjy_group/xianglilingzi/software/miniconda/bin/perl

=head1 NAME

statistics.pl - sequence information statistics

=head1 SYNOPSIS

USAGE: 

  statistics2.0.pl <input_file[,input_file1,...]>

Optional arguments:

    --input_file,-i                     Path to a sequence file
    --reflength,-r                      The length of reference sequence
    --type,-t                           The type of the input file (fastq|fasta, default=auto)
    --together                          Count all the input files as a whole
    --head                              Export header
    --help,-h                           Help message

=head1 OPTIONS

B<--input_file,-i>

    A sequence file.

B<--reflength,-r>

    The length of the reference sequence.

B<--type,-t>

    The type of the input file. The default value is auto.

B<--together>

    If you specify this option, all of the input files will be considered as a whole.

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

my (%options, @input, $name, %length, %GC, %N, @len, $N50, $N90);
my ($sumlength, $sumGC, $sumN, $seqnum, $GCp, $Np, $coverage)=(0,0,0,0,0,0,0);
GetOptions(\%options,
	'input_file|i=s'=>\@input,
	'reflength|r=i',
	'type|t=s',
	'together',
	'head',
	'help|h'
) || pod2usage();

## display documentation;
if ($options{help}){
	pod2usage({-exitval => 0, -verbose => 2, -output => \*STDERR});
}

## make sure everything passed was peachy;
&check_parameters();

foreach (@input){
	$name = &read_file($_);
	&statistics() unless (defined $options{together});
	&output($name) unless (defined $options{together});
	(%length, %GC, %N, @len) = () unless (defined $options{together});
	($sumlength, $sumGC, $sumN, $seqnum, $GCp, $Np, $coverage, $N50, $N90)= qw(0,0,0,0,0,0,0,0,0) unless (defined $options{together});
}
&statistics() if (defined $options{together});
$name =~ s/_\d\.clean\.fq\.gz//g if (defined $options{together});
&output($name) if (defined $options{together});

#########################
#######Sub Function######
#########################

sub check_parameters {
	die "you must input a --input_file" unless (@input);
}

sub read_file {
	my $file = shift;
	my $basename = basename($file);
	if ($file =~ /\.gz$/){
		open IN, "<:gzip", "$file" or die $!;
	}else{
		open IN, "$file" or die $!;
	}
	if (!defined $options{type}){
		if ($_ =~ /\.fq|\.fastq/){
			$options{type}='fastq';
		}
		if ($_ =~ /\.fa|\.fna|\.fasta/){
			$options{type}='fasta';
		}
	}
	&read_fastq if ($options{type} =~ /fastq/);
	&read_fasta if ($options{type} =~ /fasta/);
	return $basename;
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

sub statistics {
	foreach my $k (keys %length){
		$sumlength += $length{$k};
		push @len, $length{$k};
		$sumGC += $GC{$k};
		$sumN += $N{$k};
	}
	$seqnum = @len;
	$GCp = $sumGC/($sumlength-$sumN)*100;
	$Np = $sumN/$sumlength*100;
	if ($options{type} =~ /fasta/){
		my @sortlen = sort {$b<=>$a} @len;
		foreach (@sortlen){
			my $sum += $_;
			$N50 = $_ if (!$N50 && $sum > $sumlength*0.5);
			$N90 = $_ if (!$N90 && $sum > $sumlength*0.9);
		}
	}
	$coverage = $sumlength/$options{reflength} if ($options{type} =~ /fastq/ && defined $options{reflength});
}

sub output{
	my $filename = shift;
	if ($options{type} =~ /fasta/){
		printf STDOUT "%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n", 'File name', 'Seq number', 'Total length', 'GC content(%)', 'Gap rate(%)', 'N50 length', 'N90 length' if (defined $options{head});
		printf STDOUT "%15s\t%15d\t%15d\t%15.2f\t%15.2f\t%15d\t%15d\n", $filename, $seqnum, $sumlength, $GCp, $Np, $N50, $N90;
	}
	if ($options{type} =~ /fastq/){
		if (defined $options{reflength}){
			printf STDOUT "%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n", 'File name', 'Read number', 'Clean Base', 'GC content(%)', 'Gap rate(%)', 'Coverage' if (defined $options{head});
			printf STDOUT "%15s\t%15d\t%15d\t%15.2f\t%15.2f\t%15.2f\n", $filename, $seqnum, $sumlength, $GCp, $Np, $coverage;
		}else{
			printf STDOUT "%15s\t%15s\t%15s\t%15s\t%15s\n", 'File name', 'Read number', 'Clean Base', 'GC content(%)', 'Gap rate(%)' if (defined $options{head});
			printf STDOUT "%15s\t%15d\t%15d\t%15.2f\t%15.2f\n", $filename, $seqnum, $sumlength, $GCp, $Np;
		}
	}
}
