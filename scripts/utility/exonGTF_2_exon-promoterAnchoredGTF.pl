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


open (GTF,">gtf$$") or die "cannot open gtf$$";
while (<STDIN>){
  chomp $_;
  unless ($_=~/$tag/){system "rm gtf$$";die "ERROR! Line:\n$_\ndo not contain the tag $tag\nPlease set the proper -tag parameter\n";};
  $_=~s/$tag/transcript_id/ unless ($tag eq 'transcript_id');
  $_ .= " gene_id \"FAKE\";" if ($_!~/gene_id\s+\"[^\"]+\"/);
  print GTF "$_\n";
}
close GTF;

my %infoGtf       = readingGTFexons("gtf$$");
my %infoPromoters = gtfExons2gtfPromoters(%infoGtf);


my %allTags;
open (F,"<gtf$$") or  die "cannot open gtf$$";
while (my $l = <F>){
  chomp $l;
  if ($l=~/transcript_id\s+\"([^\"]+)\"/){
    my $id = $1;
    unless (defined $allTags{$id}){
      my $outLine;
      $outLine .= $infoPromoters{$id}->{'chr'}     . "\t";
      $outLine .= $infoPromoters{$id}->{'source'}  . "\t";
      $outLine .= "exon"                           . "\t";
      if ($infoPromoters{$id}->{'start'} > 0) {
      	$outLine .= $infoPromoters{$id}->{'start'} . "\t";
      }
      else {
	$outLine .= 1 . "\t";
	system "echo $id >> warnings_promotersStartingBeforePosition0";
      }	
      $outLine .= $infoPromoters{$id}->{'end'}     . "\t";
      $outLine .= $infoPromoters{$id}->{'score'}   . "\t";
      $outLine .= $infoPromoters{$id}->{'strand'}  . "\t";
      $outLine .= $infoPromoters{$id}->{'frame'}   . "\t";
      $outLine .= $infoPromoters{$id}->{'group'};
      print "$outLine\n";
    }
    $allTags{$id} = 1;
    print  "$l\n";
  }

}
close F;
system "rm gtf$$";

sub gtfExons2gtfPromoters {
    my (%infoGtf) = @_;
    my %infoPromoters;
    foreach my $txID (keys %infoGtf){
      @{$infoGtf{$txID}} = sort {$a->{"start"}  <=>  $b->{"start"}} @{$infoGtf{$txID}};
      my $spatiation = 0;
      my $span       = 500;
      my ($beginning , $promoterStart , $promoterEnd);
      if ($infoGtf{$txID}[0]{'strand'} eq "+"){
	$beginning = $infoGtf{$txID}[0]{'start'};
	$promoterStart = ($beginning - $spatiation) - $span;
	$promoterEnd   = ($promoterStart + $span) -1;
      }
      elsif ($infoGtf{$txID}[0]{'strand'} eq "-"){
	my $lastExon = scalar(@{$infoGtf{$txID}}) -1;
	$beginning = $infoGtf{$txID}[$lastExon]{'end'};
	$promoterStart = ($beginning + $spatiation) +1;
	$promoterEnd   = ($promoterStart + $span) -1;
      }
      else {die "Error! $infoGtf{$txID}[0]{'strand'} strand is not accepted. Choose either + or -\n";}

      my %hash = (
		  "chr"           => $infoGtf{$txID}[0]{'chr'},
		  "source"        => $infoGtf{$txID}[0]{'source'},
		  "feature"       => "promoter",
		  "start"         => $promoterStart ,
		  "end"           => $promoterEnd ,
		  "score"         => $infoGtf{$txID}[0]{'score'},
		  "strand"        => $infoGtf{$txID}[0]{'strand'},
		  "frame"         => $infoGtf{$txID}[0]{'frame'},
		  "group"         => $infoGtf{$txID}[0]{'group'},
		  "gene_id"       => $infoGtf{$txID}[0]{'group'}
		 );
      $infoPromoters{$txID} = \%hash;
    }
    return %infoPromoters;
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



=head1 NAME

exonGTF_2_exon-promoterAnchoredGTF.pl - given an exon GTF STDIN it returns another GTF including both the exons and the promoters of the transcripts, where the promoters are anyway labeled in the third field as exons. 
You can use this script to obtain coordinates of mature transcripts with the promoter upstream that, as suggested by Eugene Koonin, can be used to enhance the lncRNA homology discovery by Blast.

=head1 SYNOPSIS

perl exonGTF_2_exon-promoterAnchoredGTF.pl  -tag name [-help]

=head1 DESCRIPTION

=item * given GTF STDIN it considers just lines having "exon" has third field

=item * then it considers the transcripts relying on the transcript_id field. This is the default. If there is any other tag to consider instead like "hitName", use the parameter -tag to specify it.

=item * in any case the output will show the transcript_id field (instead of any specified tag) and it will also include a fake gene_id if this is not found in the input.

=item * ##WARNING## The default span of the promoter is 500kb befor the starting exon (according with the strand), and the spatiation is zero, that is it is stuck to the start exon

=item * ##WARNING## Exons and promoters will come together in the output file and there won't be any distinction among them (both labeled with a "exon" label in the third field )

=head1 OPTIONS

=item B<-tag> I<tagName>

default[transcript_id]. Choose the tag to use to indicate which exons belong to what transcripts. Typically you can have tag called like "hitName"

Help

=cut
