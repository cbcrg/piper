#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $tag;
my $help                   = 0;
# parsing parameters
my $result = GetOptions ("tag=s"      => \$tag,
                         "help"       => \$help);
$help = 1 unless ($result);
pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1) if ($help);

$tag = 'transcript_id' unless (defined $tag);


my (@sortedTags , %allTags);
open (GTF,">gtf$$") or die "cannot open gtf$$";
while (<STDIN>){
  chomp $_;
  unless ($_=~/$tag/){system "rm gtf$$";die "ERROR! Line:\n$_\ndo not contain the tag $tag\nPlease set the proper -tag parameter\n";};
  $_=~s/$tag/transcript_id/ unless ($tag eq 'transcript_id');
  $_ .= " gene_id \"FAKE\";" if ($_!~/gene_id\s+\"[^\"]+\"/);

  if ($_=~/transcript_id\s+\"([^\"]+)\"/){
    my $id = $1;
    push (@sortedTags , "$id") unless (defined $allTags{$id});
    $allTags{$id} = 1;
    print GTF "$_\n";
  }
}
close GTF;


my %infoGtf       = readingGTFexons("gtf$$");
system "rm gtf$$";

my %infoIntrons = gtfExons2gtfIntrons(%infoGtf);
foreach my $id (@sortedTags){
  foreach my $intron (0..$#{$infoIntrons{$id}}){
    my $outLine;
    $outLine .= $infoIntrons{$id}[$intron]->{'chr'}     . "\t";
    $outLine .= $infoIntrons{$id}[$intron]->{'source'}  . "\t";
    $outLine .= $infoIntrons{$id}[$intron]->{'feature'} . "\t";
    $outLine .= $infoIntrons{$id}[$intron]->{'start'}   . "\t";
    $outLine .= $infoIntrons{$id}[$intron]->{'end'}     . "\t";
    $outLine .= $infoIntrons{$id}[$intron]->{'score'}   . "\t";
    $outLine .= $infoIntrons{$id}[$intron]->{'strand'}  . "\t";
    $outLine .= $infoIntrons{$id}[$intron]->{'frame'}   . "\t";
    $outLine .= $infoIntrons{$id}[$intron]->{'group'};

    print "$outLine\n";
  }
}

sub readingGTFexons {
  my ($GTFfile) = @_;
  open (IN,"<$GTFfile") or die "cannot open the GTFfile $!";
  my %infoGtf;
  foreach my $line (<IN>){
	chomp $line;
	my ($chr , $group , $frame , $strand , $score , $end , $start , $feature , $source , $gene_id , $transcript_id);
	if ($line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.*)$/){
	    $chr     = $1;
	    $source  = $2;
	    $feature = $3;
	    $start   = min ($4, $5);
	    $end     = max ($4, $5);
	    $score   = $6;
	    $strand  = $7;
	    $frame   = $8;
	    $group   = $9;
	    if ($group =~/gene_id \"([^\"]+)\"/){
		$gene_id = $1;
	    }
	    else {die "gtf file in wrong format. Impossible to find the gene_id in line:\n$line\n when it is mandatory for this format\n";}
	    if ($group =~/transcript_id \"([^\"]+)\"/){
		$transcript_id = $1;
	    }
	    else {die "gtf file in wrong format. Impossible to find the transcript_id in line:\n$line\n when it is mandatory for this format\n";}

	    next unless ($feature eq 'exon');

	    my %hash = (
			"chr"           => $chr,
			"source"        => $source,
			"feature"       => $feature,
			"start"         => $start,
			"end"           => $end,
			"score"         => $score,
			"strand"        => $strand,
			"frame"         => $frame,
			"group"         => $group,
			"gene_id"       => $gene_id,
		       );
	    push (@{$infoGtf{$transcript_id}}, \%hash);
	}
	else {die "cannot read the line $line\n";}
      }
  close IN;
  return %infoGtf;
}


sub gtfExons2gtfIntrons {
    my (%infoGtf) = @_;
    my %infoIntrons;
    foreach my $txID (keys %infoGtf){
      @{$infoGtf{$txID}} = sort {$a->{"start"}  <=>  $b->{"start"}} @{$infoGtf{$txID}};
      foreach my $exon (0..$#{$infoGtf{$txID}}){
	next if (! defined $infoGtf{$txID}[1+$exon]);
	my $ex_start = $infoGtf{$txID}[$exon]{'start'};
	my $ex_end   = $infoGtf{$txID}[$exon]{'end'};
	my $next_ex_start = $infoGtf{$txID}[1+$exon]{'start'};

	my %hash = (
			"chr"           => $infoGtf{$txID}[$exon]{'chr'},
			"source"        => $infoGtf{$txID}[$exon]{'source'},
			"feature"       => "intron",
			"start"         => $ex_end +1 ,
			"end"           => $next_ex_start -1,
			"score"         => $infoGtf{$txID}[$exon]{'score'},
			"strand"        => $infoGtf{$txID}[$exon]{'strand'},
			"frame"         => $infoGtf{$txID}[$exon]{'frame'},
			"group"         => $infoGtf{$txID}[$exon]{'group'},
			"gene_id"       => $infoGtf{$txID}[$exon]{'group'}
		       );
	push (@{$infoIntrons{$txID}}, \%hash);
      }
    }
    return %infoIntrons;
}


=head1 NAME

exonGTF_2_intronGTF.pl - given an exon GTF STDIN it returns the relative introns in GTF format

=head1 SYNOPSIS

perl exonGTF_2_intronGTF.pl  -tag name [-help]

=head1 DESCRIPTION

=item * given GTF STDIN it considers just lines having "exon" has third field

=item * then it considers the transcripts relying on the transcript_id field. This is the default. If there is any other tag to consider instead like "hitName", use the parameter -tag to specify it.

=item * in any case the output will show the transcript_id field (instead of any specified tag) and it will also include a fake gene_id if this is not found in the input.

=head1 OPTIONS

=item B<-tag> I<tagName>

default[transcript_id]. Choose the tag to use to indicate which exons belong to what transcripts. Typically you can have tag called like "hitName"

Help

=cut
