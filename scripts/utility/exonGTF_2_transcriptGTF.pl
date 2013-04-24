#!/usr/bin/perl -w
use strict;
use warnings;
#exonGTF_2_transcriptGTF.pl - given an exon GTF STDIN it returns the transcripts in GTF format
#SYNOPSIS
#perl exonGTF_2_transcriptGTF.pl  -tag name [-help]
#DESCRIPTION
#given GTF STDIN it considers just lines having "exon" has third field
#then it assemble the transcripts relying on the transcript_id field. This is the default. If there is any other tag to consider instead like "hitName", use the parameter -tag to specify it.
# in any case the output will show the transcript_id field (instead of any specified tag) and it will also include a fake gene_id if this is not found in the input.
#OPTION -tag
#default[transcript_id]. Choose the tag to use to indicate which exons belong to what transcripts. Typically you can have tag called like "hitName"



my $tag = $ARGV[0];
$tag = 'transcript_id' unless (defined $tag);


my (@sortedTags , %allTags);
open (GTF,">gtf$$") or die "cannot open gtf$$";
while (<STDIN>){
  chomp $_;
  unless ($_=~/$tag/){system "rm gtf$$";die "ERROR! Line:\n$_\ndo not contain the tag $tag\nPlease set the proper -tag parameter\n";};
  $_=~s/$tag/transcript_id/ unless ($tag eq 'transcript_id');
  $_ .= " gene_id \"FAKE\";" if ($_!~/gene_id\s+\"[^\"]+\"/);

  if ($_=~/hitName\s+\"([^\"]+)\"/){
    my $id = $1;
    push (@sortedTags , "$id") unless (defined $allTags{$id});
    $allTags{$id} = 1;
    print GTF "$_\n";
  }
}
close GTF;


my %infoExons       = readingGTFexons("gtf$$");
my %infoTranscripts = gtfExons2Transcripts(%infoExons);
system "rm gtf$$";


foreach my $id (@sortedTags){
  my $outLine;
  $outLine .= $infoTranscripts{$id}->{'chr'}     . "\t";
  $outLine .= $infoTranscripts{$id}->{'source'}  . "\t";
  $outLine .= "transcript" . "\t";
  $outLine .= $infoTranscripts{$id}->{'start'}   . "\t";
  $outLine .= $infoTranscripts{$id}->{'end'}     . "\t";
  $outLine .= $infoTranscripts{$id}->{'score'}   . "\t";
  $outLine .= $infoTranscripts{$id}->{'strand'}  . "\t";
  $outLine .= $infoTranscripts{$id}->{'frame'}   . "\t";
  $outLine .= $infoTranscripts{$id}->{'group'};

  print "$outLine\n";
}



sub readingGTFexons {
  my ($GTFfile) = @_;
  open (IN,"<$GTFfile") or die "cannot open the GTFfile $!";
  my %infoGtf;
  foreach my $line (<IN>){
	chomp $line;
	my ($chr , $group , $frame , $strand , $score , $end , $start , $feature , $source , $gene_id , $transcript_id , $hit_id);
	if ($line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.*)$/){
	    $chr     = $1;
	    $source  = $2;
	    $feature = $3;
	    $start   = min2 ($4, $5);
	    $end     = max2 ($4, $5);
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
	    if ($group =~/hitName \"([^\"]+)\"/){
		$hit_id = $1;
	    }
	    else {die "gtf file in wrong format. Impossible to find the hit_id in line:\n$line\n when it is mandatory for pipeR\n";}


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
	    push (@{$infoGtf{$hit_id}}, \%hash);
	}
	else {die "cannot read the line $line\n";}
      }
  close IN;
  return %infoGtf;
}



sub gtfExons2Transcripts{
  my (%infoGtf) = @_;
  my %infoGtfTranscript;
  foreach my $rna (keys %infoGtf){
    my $minExonStart = 1000000000000000000000000000000000000000000000000000000000000000000000000000;
    my $maxExonEnd   = 0;
    my ($chr , $group , $frame , $strand , $score , $feature , $source , $gene_id);
    foreach my $index (0..$#{$infoGtf{$rna}}){
      $minExonStart = min2 ($minExonStart, $infoGtf{$rna}->[$index]->{"start"});
      $maxExonEnd   = max2 ($maxExonEnd, $infoGtf{$rna}->[$index]->{"end"});

      $chr           = $infoGtf{$rna}->[$index]->{"chr"};
      $group         = $infoGtf{$rna}->[$index]->{"group"};
      $frame         = $infoGtf{$rna}->[$index]->{"frame"};
      $strand        = $infoGtf{$rna}->[$index]->{"strand"};
      $score         = $infoGtf{$rna}->[$index]->{"score"};
      $feature       = $infoGtf{$rna}->[$index]->{"feature"};
      $source        = $infoGtf{$rna}->[$index]->{"source"};
      $gene_id       = $infoGtf{$rna}->[$index]->{"gene_id"};
    }
    $infoGtfTranscript{$rna}->{'start'}        = $minExonStart;
    $infoGtfTranscript{$rna}->{'end'}          = $maxExonEnd;
    $infoGtfTranscript{$rna}->{'chr'}          = $chr;
    $infoGtfTranscript{$rna}->{'group'}        = $group;
    $infoGtfTranscript{$rna}->{'frame'}        = $frame;
    $infoGtfTranscript{$rna}->{'strand'}       = $strand;
    $infoGtfTranscript{$rna}->{'score'}        = $score;
    $infoGtfTranscript{$rna}->{'feature'}      = $feature;
    $infoGtfTranscript{$rna}->{'source'}       = $source;
    $infoGtfTranscript{$rna}->{'gene_id'}      = $gene_id;
  }
  return %infoGtfTranscript;
}


sub min2 {
  my ($a, $b) = @_;

  return $b unless (defined $a);
  return $a unless (defined $b);
  return (($a < $b)? $a: $b);
}

sub max2 {
  my ($a, $b) = @_;

  return $b unless (defined $a);
  return $a unless (defined $b);
  return (($a > $b)? $a: $b);
}

