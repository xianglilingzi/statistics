#!/public/agis/znjy_group/xianglilingzi/software/miniconda/bin/perl

=head1 NAME

alignment_statistics1.0.pl - alignments information statistics

=head1 SYNOPSIS

USAGE: 

  alignment_statistics1.0.pl

Optional arguments:

    --unrmdup_flagstat_file,-uf                   Path to a flagstat file before deleting duplicates
    --markdup_flagstat_file,-mf                   Path to a flagstat file after deleting duplicates
    --unrmdup_depth_file,-ud                      Path to a depthstat file before deleting duplicates
    --markdup_depth_file,-md                      Path to a depthstat file after deleting duplicates
    --ref_length,-rl                              The length of a reference sequence
    --head                                        Export title
    --help,-h                                     Help message

=head1 OPTIONS

B<--unrmdup_flagstat_file,-uf>

    Path to a flagstat file before deleting secondary alignments, supplementary alignments and PCR duplicates.

B<--markdup_flagstat_file,-mf>

    Path to a flagstat file after deleting secondary alignments, supplementary alignments and PCR duplicates.

B<--unrmdup_depth_file,-ud>

    Path to a depthstat file before deleting secondary alignments, supplementary alignments and PCR duplicates.

B<--markdup_depth_file,-md>

    Path to a depthstat file after deleting secondary alignments, supplementary alignments and PCR duplicates.

B<--ref_length,-rl>

    The length of a reference sequence.

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
use File::Basename;

my (%options);
GetOptions(\%options,
	'unrmdup_flagstat_file|uf=s',
	'markdup_flagstat_file|mf=s',
	'unrmdup_depth_file|ud=s',
	'markdup_depth_file|md=s',
	'ref_length|rl=i',
	'head',
	'help|h'
) || pod2usage();

## display documentation;
if ($options{help}){
	pod2usage({-exitval => 0, -verbose => 2, -output => \*STDERR});
}

my ($uread, $umap)=&read_flagstat($options{unrmdup_flagstat_file}) if ($options{unrmdup_flagstat_file});
my ($mread, $mmap)=&read_flagstat($options{markdup_flagstat_file}) if ($options{markdup_flagstat_file});
my $umapp = $umap/$uread*100 if ($options{unrmdup_flagstat_file});
my $dupp = $mread/$uread*100 if ($options{markdup_flagstat_file} && $options{unrmdup_flagstat_file});
my $mmapp = $mmap/$mread*100 if ($options{markdup_flagstat_file} && $options{unrmdup_flagstat_file});

my ($uavdep, $ucov1p, $ucov2p, $ucov3p, $ucov4p)=&read_depth($options{unrmdup_depth_file}) if ($options{unrmdup_depth_file});
my ($mavdep, $mcov1p, $mcov2p, $mcov3p, $mcov4p)=&read_depth($options{markdup_depth_file}) if ($options{markdup_depth_file});

(my $filename = basename($options{markdup_flagstat_file})) =~ s/.markdup.flagstat//g;

printf STDOUT "%15s\t", 'File name' if($options{head});
printf STDOUT "%15s\t%15s\t%15s\t", 'Clean Read', 'C Mapping Read', 'C Mapping Rate' if ($options{unrmdup_flagstat_file} && $options{head});
printf STDOUT "%15s\t%15s\t%15s\t", 'C Average Depth', 'C Coverage1X', 'C Coverage2X', 'C Coverage3X', 'C Coverage4X' if ($options{unrmdup_depth_file} && $options{head});
printf STDOUT "%15s\t%15s\t%15s\t%15s\t", 'Rmdup Read', 'Rmdup Rate', 'R Mapping Read', 'R Mapping Rate' if ($options{unrmdup_depth_file} && $options{markdup_flagstat_file} && $options{head});
printf STDOUT "%15s\t%15s\t%15s\t", 'R Average Depth', 'R Coverage1X', 'R Coverage2X', 'R Coverage3X', 'R Coverage4X' if ($options{markdup_depth_file} && $options{head});
printf STDOUT "\n" if (defined $options{head});
printf STDOUT "%15s\t", $filename;
printf STDOUT "%15d\t%15d\t%15.2f\t", $uread, $umap, $umapp if ($options{unrmdup_flagstat_file});
printf STDOUT "%15.2f\t%15.2f\t%15.2f\t%15.2f\t%15.2f\t", $uavdep, $ucov1p, $ucov2p, $ucov3p, $ucov4p if ($options{unrmdup_depth_file});
printf STDOUT "%15d\t%15.2f\t%15d\t%15.2f\t", $mread, $dupp, $mmap, $mmapp if ($options{markdup_flagstat_file} && $options{unrmdup_depth_file});
printf STDOUT "%15.2f\t%15.2f\t%15.2f\t%15.2f\t%15.2f\t", $mavdep, $mcov1p, $mcov2p, $mcov3p, $mcov4p if ($options{markdup_depth_file});
printf STDOUT "\n";


###########################
########Sub function#######
###########################

sub read_flagstat {
	my $file = shift @_;
	my ($read, $pair, $single)=(0,0,0);
	open IN, "$file" or die $!;
	while (<IN>){
		if (/(\d+)\s+\+\s+\d+\s+read\d\s+/){
			$read+=$1;
		}
		if (/(\d+)\s+\+\s+\d+\s+with\s+itself\s+and\s+mate\s+mapped\s+/){
			$pair+=$1;
		}
		if (/(\d+)\s+\+\s+\d+\s+singletons\s+\(\d+\.\d+\%\s+:\s+N\/A\)\s+/){
			$single+=$1;
		}
	}
	close IN;
	my $map=$pair+$single;
	return ($read, $map);
}

sub read_depth {
	my $file = shift @_;
	my ($depth, $cov1, $cov2, $cov3, $cov4, $base)=(0,0,0,0,0,0);
	open IN, "$file" or die $!;
	while (<IN>){
		if (/\S+\s+\d+\s+(\d+)\s+/){
			$depth+=$1;
			$cov1++ if ($1>0);
			$cov2++ if ($1>1);
			$cov3++ if ($1>2);
			$cov4++ if ($1>3);
			$base++;
		}
	}
	close IN;
	$base = $options{ref_length} if ($options{ref_length});
	my $avdep = $depth/$base;
	my $cov1p = $cov1/$base*100;
	my $cov2p = $cov2/$base*100;
	my $cov3p = $cov3/$base*100;
	my $cov4p = $cov4/$base*100;
	return ($avdep, $cov1p, $cov2p, $cov3p, $cov4p);
}
