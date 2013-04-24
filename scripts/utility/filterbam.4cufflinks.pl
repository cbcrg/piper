#!/soft/bin/perl

#    GRAPE
#    Copyright (C) 2011 Centre for Genomic Regulation (CRG)
#
#    This file is part of GRAPE.
#
#    GRAPE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    GRAPE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GRAPE.  If not, see <http://www.gnu.org/licenses/>.

# Author : David Gonzalez, david.gonzalez@crg.eu


use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/[^\/]*$/lib/;
    unshift @INC, "$libdir";
}

# Objective
# This script should read a bam file, extract the sam, sort it and produce an
# output file contaiing only the uniquely mapping reads in it

use Getopt::Long;
use Bio::DB::Fasta;
use Bio::DB::Sam;
use Pod::Usage;
use RNAseq_pipeline3 ('get_fh','run_system_command');
use RNAseq_pipeline_settings3 ('read_config_file');
#use Tools::Bam ('generate_sorted_bam');
use Data::Dumper;

# Declare variables & Get command line options
my $infile;
my $tmpdir;
my $genomefile='';
my $reindex=0;

my %options=%{read_config_file()};
$tmpdir=$options{'LOCALDIR'};

GetOptions('genome|g=s' => \$genomefile);

pod2usage("Reference file (-g) required") unless ($genomefile);

# Build the subroutine that will get the sequence
*get_sequence=get_chr_subseq($genomefile,
			     $reindex);

# Get the infile
$infile=shift;
unless($infile) {
    die "No input file provided\n";
}

# Process the file
print STDERR "Processing $infile\n";
print STDERR "Extracting uniquely mapping reads acordig to the alignmnent\n";

# produce a sam file sorted by read_id, pipe it into the script and print out
# only those cases where the read id appears only once.
my $outfile=filterbam($infile);

# Get the header from the original bam file
my $headerfile=get_header($infile);

# Concatenate both files and create the sorted bam
my $bamfile=build_output_bam($headerfile,
			     $outfile,
			     $infile);

# Clean up
my $command="rm $headerfile $outfile";
run_system_command($command);

exit;

sub build_output_bam {
    my $header=shift;
    my $samfile=shift;
    my $bamfile=shift;

    $bamfile=~s/bam$/filtered/;
    $bamfile=~s/.*\///;

    my $command="cat $header $samfile | samtools view -S -b - |samtools sort - $bamfile";
    run_system_command($command);

    
}

sub get_header {
    my $infile=shift;

    my $outfile=$infile;

    $outfile=~s/.*\///;
    $outfile=$$.'.header.'.$outfile;

    print STDERR "Extracting header from $infile\n";

    my $command="samtools view -H $infile > $outfile";
    run_system_command($command);

    return($outfile);
}

sub filterbam {
    my $infile=shift;
    my $outfile=$infile;

    $outfile=~s/.*\///;
    $outfile=$$.'.'.$outfile;

    my $command="samtools view $infile|";
    print STDERR $command,"\n";

    my $infh;
    open($infh,$command);

    my $outfh=get_fh($outfile,1);

    my $old_read_id='';
    my @lines=();

    my %strand=("GTAG" => '+',
		"CTAC" => '-',
		"GCAG" => '+',
		"CTGC" => '-',
		"ATA." => '+',
		".TAT" => '-',
		"GTAT" => '+',
		"ATAC" => '-');

    while (my $line=<$infh>) {
	chomp($line);
	my @line=split("\t",$line);

	my $read_id=$line[0];
	my $flag=$line[1];
	my $chr=$line[2];
	my $start=$line[3];
	my $cigar=$line[5];
	my $flags=$line[11];
	my $cuff_flag='XS:A:+';

	my @splits=split('N',$cigar);
	my ($start1,$start2,$end1,$end2,$key);
	if (@splits > 2) {
	    print STDERR $cigar,"\n";
	} elsif (@splits == 2) {
	    # The read is split
	    my @coords=split('M',$splits[0]);
	    unless(@coords == 2) {
		die "Wrong number of coords\n";
	    }
	    $splits[1]=~s/M//;
	    $start1=$start + $coords[0];
	    $end1=$start + $coords[0] + 1;
	    $start2= $start + $coords[0] + $coords[1] - 2 ;
	    $end2= $start + $coords[0] + $coords[1] - 1;


	    my $seq1=get_sequence($chr,$start1,$end1);
	    my $seq2=get_sequence($chr,$start2,$end2);

	    $key=$seq1.$seq2;

	    if ($strand{$key}) {
		$cuff_flag='XS:A:'.$strand{$key};
	    }
	    if ($flags) {
		$flags=join(' ',$flags,$cuff_flag);
	    } else {
		$flags=$cuff_flag;
	    }
	    $line=join("\t",@line[0..10],$flags);
	}
	$line.="\n";

	# If the last previous read is differnt check if it is unique
	# If it is unique print it
	if ($read_id ne $old_read_id) {
	    if (@lines == 1) {
		
		print $outfh @lines;
	    }
	    $old_read_id=$read_id;
	    @lines=();
	}
#	print $line;
	push @lines,$line;
    }

    # Print the last record
    if (@lines == 1) {
	print $outfh @lines;
    }
    close($outfh);
    close($infh);

    return($outfile);
}

sub get_chr_subseq {
    my $genomefile=shift;
    my $reindex=shift;

    my %chromosomes;
    my %cache;

    if ($reindex) {
	print STDERR "Regenerating index...\n";
    }

    # Generate the chromosomes
    print STDERR "Generating $genomefile index...";
    # Here we reindex because the script seems not to find the required
    # index if not when the index is in a different dir
    my $db=Bio::DB::Fasta->new($genomefile,
			       -reindex => $reindex);
    print STDERR "done\n";

    # Build the subroutine for the search
    my $get_seq_from_chrom= sub {
	my $chr=shift;
	my $start=shift;
	my $end=shift;

	# Get the slice
	my $seq;
	$seq=$db->subseq($chr,$start,$end);
    };
    return($get_seq_from_chrom);
}

=head1 NAME
    
    filterbam.4cufflinks.pl
    
=head1 SYNOPSIS
    
filterbam.4cufflinks.pl [-h][-man][-genome <file>] file.bam
    
  Options:
    -help|h          Brief help message
    -man             Full documentation
    -genome|g file   Indicate the name of the reference file to use
        

    This script will take a bam file and add to thoe split reads the cufflinks
    flag XS:A:+ or XS:A:- This is done in the same manner as tophat does, by 
    looking at the two nucleotides at the end of the intron at both ends and
    assigning the strand accordingly. If the nucleotides do not comply the
    strand is set to + by default.

=head1 OPTIONS
    
=over 8
    
=item B<-help|h>
    
    Print a brief help message and exits.
    
=item B<-man>
    
    Prints the manual page and exits.

=item B<-debug|d>
    
    Currently only exits after retrieving the first sequence

=item B<-genome|g> file
    
    Specify the reference file from which the sequences should be retrieved. The
    format of this file should be multifasta.

=item B<-border> n
    
    Instead of retrieving the region between the coordinates in the gff entry
    return the n flanking sequences on both sides

=item B<-reindex|r>
    
    Regenerate the Bio::Fasta::DB index for fast sequence retrieval, this
    should not be necessary unless the index cannot be read correctly

=item B<-minus>
    
    Return the sequence on the minus strand for those gff regions which contain
    no strand information (.), by default the plus sequence is returned

=back
    
=head1 DESCRIPTION
    
    This program should take a gff file and for each of the records extract the
    sequence from a reference sequence file

=cut
