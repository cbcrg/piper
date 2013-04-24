#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use File::Basename;
use Switch;
#Author: Giovanni Bussotti

#HELP
my $help  = 0;
foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '-help') or ($ARGV[$field] eq '--help'));
}
if ($help > 0){
 help_message();
}

#TAKE OPTIONS
acceptedVariableSpace();
my ($type , $repeatThreshold , $chr_subseq , $promoters_size , $feature2extract  ,  $experimentName , $pipelineDirName) = options();
my $results_out = "${pipelineDirName}/experiments/${experimentName}/results/";
my $outDir      = "${pipelineDirName}/experiments/${experimentName}/qualityCheck/extractedFeatues/";

#EXONS
if (($feature2extract eq 'exons') or ($feature2extract eq 'all')){
  if (-d "${outDir}/exons_${type}_rep$repeatThreshold"){(system "rm -rf ${outDir}/exons_${type}_rep$repeatThreshold") == 0 or die "Error[extractFeatures.pl]! Cannot remove ${outDir}/exons_${type}_rep$repeatThreshold\n";}
  (system "mkdir -p ${outDir}/exons_${type}_rep$repeatThreshold") == 0 or die "Error[extractFeatures.pl]! cannot execute\nmkdir -p ${outDir}/exons_${type}_rep$repeatThreshold\n$!\n";
  splitySpeciesByQuery("exons_${type}_rep$repeatThreshold");
  GTF_TO_FASTA_transcript("${outDir}/exons_${type}_rep$repeatThreshold/","no");
}

#PROMOTERS
if (($feature2extract eq 'promoters') or ($feature2extract eq 'all')){
  if (-d "${outDir}/promoters_${type}_rep$repeatThreshold"){(system "rm -rf ${outDir}/promoters_${type}_rep$repeatThreshold") == 0 or die "Error[extractFeatures.pl]! Cannot remove ${outDir}/promoters_${type}_rep$repeatThreshold\n";}
  (system "mkdir -p ${outDir}/promoters_${type}_rep$repeatThreshold") == 0 or die "Error[extractFeatures.pl]! cannot execute\nmkdir -p ${outDir}/promoters_${type}_rep$repeatThreshold\n$!\n";
  #TAKE GTF
  opendir (EO,"$results_out") or die "Error[extractFeatures.pl]! Cannot open $results_out\n$!\n";
  while (my $exGtfFile = readdir(EO)){
    next if (($exGtfFile eq '.') || ($exGtfFile eq '..') || ($exGtfFile=~/\.fa$/) || ($exGtfFile !~/\.${type}\.rep${repeatThreshold}\.ex\.gtf$/));
    my @sortedTags    = takeSortedTags("${results_out}/${exGtfFile}");
    my %infoGtf       = readingGTFexons("${results_out}/${exGtfFile}");
    my %infoPromoters = gtfExons2gtfPromoters(%infoGtf);
    printGtf(\@sortedTags,\%infoPromoters,'promoters');
  }
  closedir EO;
  #TAKE FASTA
  GTF_TO_FASTA_transcript("${outDir}/promoters_${type}_rep$repeatThreshold/","no");
}

#EXON GTF TO INTRONS
if (($feature2extract eq 'introns') or ($feature2extract eq 'all')){
  if (-d "${outDir}/introns_${type}_rep$repeatThreshold"){(system "rm -rf ${outDir}/introns_${type}_rep$repeatThreshold") == 0 or die "Error[extractFeatures.pl]! Cannot remove ${outDir}/introns_${type}_rep$repeatThreshold\n";}
  (system "mkdir -p ${outDir}/introns_${type}_rep$repeatThreshold") == 0 or die "Error[extractFeatures.pl]! cannot execute\nmkdir -p ${outDir}/introns_${type}_rep$repeatThreshold\n$!\n";
  #TAKE GTF
  opendir (EO,"$results_out") or die "Error[extractFeatures.pl]! Cannot open $results_out\n$!\n";
  while (my $exGtfFile = readdir(EO)){
    next if (($exGtfFile eq '.') || ($exGtfFile eq '..') || ($exGtfFile=~/\.fa$/) || ($exGtfFile !~/\.${type}\.rep${repeatThreshold}\.ex\.gtf$/));
    my @sortedTags    = takeSortedTags("${results_out}/${exGtfFile}");
    my %infoGtf       = readingGTFexons("${results_out}/${exGtfFile}");
    my %infoIntrons   = gtfExons2gtfIntrons(%infoGtf);
    printGtf(\@sortedTags,\%infoIntrons,'introns');
  }
  closedir EO;
  #TAKE FASTA
  GTF_TO_FASTA_transcript("${outDir}/introns_${type}_rep$repeatThreshold/","yes");
}

#EXON GTF TO TRANSCRIPTS
if (($feature2extract eq 'transcripts') or ($feature2extract eq 'all')){
  if (-d "${outDir}/transcripts_${type}_rep$repeatThreshold"){(system "rm -rf ${outDir}/transcripts_${type}_rep$repeatThreshold") == 0 or die "Error[extractFeatures.pl]! Cannot remove ${outDir}/transcripts_${type}_rep$repeatThreshold\n";}
  (system "mkdir -p ${outDir}/transcripts_${type}_rep$repeatThreshold") == 0 or die "Error[extractFeatures.pl]! cannot execute\nmkdir -p ${outDir}/transcripts_${type}_rep$repeatThreshold\n$!\n";
  #TAKE FASTA   you have to do it first for the special case of transcript. This way it is correct
  splitySpeciesByQuery("transcripts_${type}_rep$repeatThreshold");
  GTF_TO_FASTA_transcript("${outDir}/transcripts_${type}_rep$repeatThreshold/","yes");
  (system "rm -rf ${outDir}/transcripts_${type}_rep$repeatThreshold/*gtf") == 0 or die "Error[extractFeatures.pl]! Cannot execute:\nrm -rf ${outDir}/transcripts_${type}_rep$repeatThreshold/*gtf  \n$!\n";
  #TAKE GTF
  opendir (EO,"$results_out") or die "Error[extractFeatures.pl]! Cannot open $results_out\n$!\n";
  while (my $exGtfFile = readdir(EO)){
    next if (($exGtfFile eq '.') || ($exGtfFile eq '..') || ($exGtfFile=~/\.fa$/));
    my @sortedTags    = takeSortedTags("${results_out}/${exGtfFile}");
    my %infoGtf       = readingGTFexons("${results_out}/${exGtfFile}");
    my %infoTranscripts = gtfExons2Transcripts(%infoGtf);
    printGtf(\@sortedTags,\%infoTranscripts,'transcripts');
  }
  closedir EO;
}












#FUNCTIONS
sub options {
  my ($type , $repeatThreshold , $chr_subseq , $promoters_size , $feature2extract , $experimentName , $pipelineDirName);
  my $spyExperiment           = 1;
  my $spyPipelineDir          = 1;
  my $spyPromoters_size       = 1;
  my $spyFeature2extract      = 1;
  my $spyChr_subseq           = 1;
  my $spyRepeatThreshold      = 1;
  my $spyType                 = 1;

  foreach my $field (0..$#ARGV){
    if ($ARGV[$field] eq '-type'){
	$type  = $ARGV[1+$field];
	$spyType = 2;
	next;
    }
    if ($spyType == 2){
	$spyType = 3;
	next;
    }
    if ($ARGV[$field] eq '-repeatThreshold'){
	$repeatThreshold  = $ARGV[1+$field];
	$spyRepeatThreshold = 2;
	next;
    }
    if ($spyRepeatThreshold == 2){
	$spyRepeatThreshold = 3;
	next;
    }
    if ($ARGV[$field] eq '-chr_subseq'){
	$chr_subseq = $ARGV[1+$field];
	$spyChr_subseq = 2;
	next;
    }
    if ($spyChr_subseq == 2){
	$spyChr_subseq = 3;
	next;
    }
    if ($ARGV[$field] eq '-feature'){
	$feature2extract = $ARGV[1+$field];
	$spyFeature2extract = 2;
	next;
    }
    if ($spyFeature2extract == 2){
	$spyFeature2extract = 3;
	next;
    }
    if ($ARGV[$field] eq '-promoters_size'){
	$promoters_size = $ARGV[1+$field];
	$spyPromoters_size = 2;
	next;
    }
    if ($spyPromoters_size == 2){
	$spyPromoters_size = 3;
	next;
    }
    if ($ARGV[$field] eq '-experiment'){
	$experimentName = $ARGV[1+$field];
	$spyExperiment = 2;
	next;
    }
    if ($spyExperiment == 2){
	$spyExperiment = 3;
	next;
    }
    if ($ARGV[$field] eq '-pipeline_dir'){
      $pipelineDirName = $ARGV[1+$field];
      $spyPipelineDir = 2;
      next;
    }
    if ($spyPipelineDir == 2){
      $spyPipelineDir = 3;
      next;
    }

  }
  die "Error[extractFeature.pl]! Please provide the -experiment parameter\n"      if (! defined $experimentName);
  die "Error[extractFeature.pl]! Please provide the -pipeline_dir parameter\n"    if (! defined $pipelineDirName);
  die "Error[extractFeature.pl]! Please provide the -repeatThreshold parameter\n" if (! defined $repeatThreshold);
  die "Error[extractFeature.pl]! Please provide the -type parameter\n"            if (! defined $type);

  #defaults
  $feature2extract = 'all'           if (! defined $feature2extract);
  $promoters_size  = 500             if (! defined $promoters_size);
  $chr_subseq      = 'chr_subseq'    if (! defined $chr_subseq);
  return ($type , $repeatThreshold , $chr_subseq , $promoters_size , $feature2extract ,  $experimentName , $pipelineDirName);
}

sub acceptedVariableSpace {
  my %space = ('-type' => 1 , '-repeatThreshold' => 1 , '-experiment' => 1 , '-pipeline_dir' => 1   , '-feature' => 1 , '-promoters_size' => 1 , '-chr_subseq' => 1);
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[extractFeature.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}

sub help_message {
my $helpMessage = "\nNAME
extractFeatures.pl - it extracts specific features of the results files. The features supported so far are promoters, introns and transcripts

SYNOPSIS
extractFeatures.pl -experiment -pipeline_dir -repeatThreshold -type [-promoters_size -chr_subseq -feature]

DESCRIPTION
   * By default it extracts all the features
   * By default the promoters are considered the region 500 nucleotide upstream the transcription starting site (TSS)
   * By default the intronic sequences are all concatenated in one single fasta for each transcript (but you can change this by edit the printFunc function)
CAVEATS:
Promoters might be truncated if their start position is before the chromosome start. In this case the promoter start will be arranged to 1 and a warning file will be generated.
If a promoter passes the chromosome end, the gtf file will not be amended, but the FASTA of the promoter will be anyway truncated (crh_subseq stops where the chr stops. No error message is returned, but just warnings)

";
print "$helpMessage\n\n\n";
exit;
}


sub readingGTFexons {
  my ($GTFfile) = @_;
  open (IN,"<$GTFfile") or die "Error[extractFeatures.pl]! cannot open the GTFfile $!";
  my %infoGtf;
  my $species;
  my $baseName = basename($GTFfile);
  if ($baseName=~/(\S+)\.${type}\.rep${repeatThreshold}\.ex\.gtf$/){
    $species = $1;
  }
  else {die "Error[extractFeatures.pl]! Cannot take the species name";}
  foreach my $line (<IN>){
	chomp $line;
	my ($chr , $group , $frame , $strand , $score , $end , $start , $feature , $source , $gene_id , $hitName);
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
	    $group .= " ___species___ \"$species\";";
	    if ($group =~/gene_id \"([^\"]+)\"/){
		$gene_id = $1;
	    }
	    else {die "Error[extractFeatures.pl]! gtf file in wrong format. Impossible to find the gene_id in line:\n$line\n when it is mandatory for this format\n";}
	    if ($group =~/hitName \"([^\"]+)\"/){
		$hitName = $1;
	    }
	    else {die "Error[extractFeatures.pl]! gtf file in wrong format. Impossible to find the hitName in line:\n$line\n when it is mandatory for this format\n";}

	    next unless ($feature eq 'exon');
	    $group =~s/"/\\"/g;

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
	    push (@{$infoGtf{$hitName}}, \%hash);
	}
	else {die "Error[extractFeatures.pl]! cannot read the line $line\n";}
      }
  close IN;
  return %infoGtf;
}


sub gtfExons2gtfPromoters {
    my (%infoGtf) = @_;
    my %infoPromoters;
    foreach my $txID (keys %infoGtf){
      @{$infoGtf{$txID}} = sort {$a->{"start"}  <=>  $b->{"start"}} @{$infoGtf{$txID}};
      my $spatiation = 0;
      my $span       = $promoters_size;
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
      else {die "Error[extractFeatures.pl]! $infoGtf{$txID}[0]{'strand'} strand is not accepted. Choose either + or -\n";}

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


sub getReverseComplement {
  my ($sequence) = @_;
  $sequence = reverse $sequence;
  my $complementedSequence = "";
  for (my $key = 0; $key < length($sequence); $key++) {
    my $char = substr $sequence, $key, 1;
    switch ($char) {
      case "A" {
        $complementedSequence .= "T";
      }
      case "C" {
        $complementedSequence .= "G";
      }
      case "G" {
        $complementedSequence .= "C";
      }
      case "T" {
        $complementedSequence .= "A";
      }
      case "U" {
        $complementedSequence .= "A";
      }
      case "a" {
        $complementedSequence .= "t";
      }
      case "c" {
        $complementedSequence .= "g";
      }
      case "g" {
        $complementedSequence .= "c";
      }
      case "t" {
        $complementedSequence .= "a";
      }
      case "u" {
        $complementedSequence .= "a";
      }
      case "(" {
        $complementedSequence .= ")";
      }
      case ")" {
        $complementedSequence .= "(";
      }
      case "<" {
        $complementedSequence .= ">";
      }
      case ">" {
        $complementedSequence .= "<";
      }
      else {
        $complementedSequence .= $char;
      }
    }
  }
  return $complementedSequence;
}




sub string2FASTA {
    my ($string) = @_;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    my @fake = split (//, "$string");
    my $count = 0;
    my @copy;
    foreach my $letter (@fake){
        $count++;
        push (@copy,$letter);
        if ($count == 100){
            push (@copy, "\n");
            $count = 0;
        }
    }
    chomp $string if ($count == 0);
    $string       = join( "", @copy );
    return ($string);
}


sub takeSortedTags {
  my ($inputGtf) = @_;
  open (E, "<$inputGtf") or die "Error[extractFeatures.pl]! Cannot open $inputGtf\n$!\n";
  my (@sortedTags , %allTags);
  while (<E>){
    if ($_=~/hitName\s+\"([^\"]+)\"/){
      my $id = $1;
      push (@sortedTags , "$id") unless (defined $allTags{$id});
      $allTags{$id} = 1;
    }
    else{die "Error[extractFeatures.pl]! line\n$_\nin file $inputGtf does not contain the hitName field\n";}
  }
  close E;
  return @sortedTags;
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

sub gtfExons2Transcripts{
  my (%infoGtf) = @_;
  my %infoGtfTranscript;
  foreach my $rna (keys %infoGtf){
    my $minExonStart = 1000000000000000000000000000000000000000000000000000000000000000000000000000;
    my $maxExonEnd   = 0;
    my ($chr , $group , $frame , $strand , $score , $feature , $source , $gene_id);
    foreach my $index (0..$#{$infoGtf{$rna}}){
      $minExonStart = min ($minExonStart, $infoGtf{$rna}->[$index]->{"start"});
      $maxExonEnd   = max ($maxExonEnd, $infoGtf{$rna}->[$index]->{"end"});

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


sub splitySpeciesByQuery {
  my ($outFeatureFolder) = @_;
  opendir (EO,"$results_out") or die "Error[extractFeatures.pl]! Cannot open $results_out\n$!\n";
  while (my $exGtfFile = readdir(EO)){
    next if (($exGtfFile eq '.') || ($exGtfFile eq '..') || ($exGtfFile=~/\.fa$/) || ($exGtfFile !~/\.${type}\.rep${repeatThreshold}\.ex\.gtf$/));
    my $species;
    if ($exGtfFile=~/(\S+)\.${type}\.rep${repeatThreshold}\.ex\.gtf$/){
      $species = $1;
    }
    else {die "Error[extractFeatures.pl]! Cannot take the species name\n$!\n";}
    open (E, "<${results_out}/${exGtfFile}") or die "Error[extractFeatures.pl]! Cannot open ${results_out}/${exGtfFile}\n$!\n";
    while (<E>){
      chomp $_;
      $_ .= " ___species___ \"$species\";";
      if ($_=~/hitName\s+\"([^\"]+)\"/){
	my $id = $1;
	if ($_=~/^(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+)(.*)$/){
	  my $a = $1;
	  my $b = $2;
	  $b =~s/"/\\"/g;
	  $_ = $a . $b;
	}
	(system "echo \"$_\" >> ${outDir}/${outFeatureFolder}/${id}.ex.gtf") == 0 or die "Error[extractFeatures.pl]! Cannot execute\necho \"$_\" >> ${outDir}/${outFeatureFolder}/${id}.ex.gtf\n$!\n";
      }
      else{die "Error[extractFeatures.pl]! line\n$_\nin file $exGtfFile does not contain the hitName field\n";}
    }
    close E;
  }
  closedir EO;
}


sub printGtf {
  my ($sortedTags_p , $info_p  , $kind) = @_;
  my @sortedTags = @{$sortedTags_p};
  my %info       = %{$info_p};
  my ($f,$t);
  if ($kind eq 'promoters')  {$f = "promoters_${type}_rep$repeatThreshold"  ; $t = 'pr';}
  if ($kind eq 'transcripts'){$f = "transcripts_${type}_rep$repeatThreshold"; $t = 'tx';}
  if ($kind eq 'introns')    {$f = "introns_${type}_rep$repeatThreshold"    ; $t = 'in';}

  if (($kind eq 'promoters') or ($kind eq 'transcripts')){
    foreach my $id (@sortedTags){
      my $outLine;
      $outLine .= $info{$id}->{'chr'}     . "\t";
      $outLine .= $info{$id}->{'source'}  . "\t";
      $outLine .= $info{$id}->{'feature'} . "\t";
      if ($info{$id}->{'start'} > 0) {
      	$outLine .= $info{$id}->{'start'} . "\t";
      }
      else {
	$outLine .= 1 . "\t";
	system "echo $id >> warnings_promotersStartingBeforePosition0";
      }
      $outLine .= $info{$id}->{'end'}     . "\t";
      $outLine .= $info{$id}->{'score'}   . "\t";
      $outLine .= $info{$id}->{'strand'}  . "\t";
      $outLine .= $info{$id}->{'frame'}   . "\t";
      $outLine .= $info{$id}->{'group'}   . "\t";
      system "echo \"$outLine\" >> ${outDir}/${f}/${id}.${t}.gtf"
    }
  }
  if ($kind eq 'introns'){
    foreach my $id (@sortedTags){
      foreach my $intron (0..$#{$info{$id}}){
	my $outLine;
	$outLine .= $info{$id}[$intron]->{'chr'}     . "\t";
	$outLine .= $info{$id}[$intron]->{'source'}  . "\t";
	$outLine .= $info{$id}[$intron]->{'feature'} . "\t";
	$outLine .= $info{$id}[$intron]->{'start'}   . "\t";
	$outLine .= $info{$id}[$intron]->{'end'}     . "\t";
	$outLine .= $info{$id}[$intron]->{'score'}   . "\t";
	$outLine .= $info{$id}[$intron]->{'strand'}  . "\t";
	$outLine .= $info{$id}[$intron]->{'frame'}   . "\t";
	$outLine .= $info{$id}[$intron]->{'group'}   . "\t";
	system "echo \"$outLine\" >> ${outDir}/${f}/${id}.${t}.gtf"
      }
    }
  }
}


sub GTF_TO_FASTA_transcript{
  my ($gtfDir , $concatenate) = @_;
  my (%transcriptInfo);
  opendir (G,"$gtfDir") or die "Error[extractFeatures.pl]! Cannot open $gtfDir\n$!\n";
  while (my $gtfFile = readdir(G)){
    next if (($gtfFile eq '.') || ($gtfFile eq '..') || ($gtfFile=~/\.fa$/));
    my $outFastaFile = $gtfFile;
    $outFastaFile =~s/\.\S+\.\S+$/\.fa/;
    open (GTF,"<${gtfDir}/$gtfFile") or die "Error[extractFeatures.pl]!cannot open ${gtfDir}/${gtfFile} :$!\n";
    foreach my $line (<GTF>){
      my ($chr , $start , $end , $strand , $geneId , $hitName);
      if ($line =~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+/){
	$chr          = $1;
	$start        = min ($2, $3);
	$end          = max ($3, $2);
	$strand       = $4;
      }
      else {die "Error[extractFeatures.pl]! cannot parse the line $line \n";}
      if ($line=~/gene_id\s+"([^"]+)"/){
	$geneId       = $1;
      }
      else {die "Error[extractFeatures.pl]! cannot parse the line $line \n";}
      if ($line=~/hitName\s+"([^"]+)/){
	$hitName = $1;
      }
      my $species;
      if ($line=~/___species___ \"([^\"]+)\";/){
	$species = $1;
      }
      $transcriptInfo{$hitName}{$species}{'chr'}    = $chr;
      $transcriptInfo{$hitName}{$species}{'strand'} = $strand;
      my %exon = ("start" => $start , "end" => $end);
      push (@{$transcriptInfo{$hitName}{$species}{'exons'}} , \%exon);
    }
    close GTF;
  }


  #taking the FASTA of each exon
  printFunc(\%transcriptInfo,$gtfDir,$concatenate);

}



sub printFunc {
  my ($transcriptInfo_p , $gtfDir , $concatenate) = @_;
  my %transcriptInfo = %{$transcriptInfo_p};

  foreach my $hitName (keys %transcriptInfo){
    open (OUT,">${gtfDir}/${hitName}.fa") or die "Error[extractFeatures.pl]!cannot open ${gtfDir}/${hitName}.fa :$!\n";
    foreach my $species (keys %{$transcriptInfo{$hitName}}){
      my $fileChr = "${pipelineDirName}/allGenomeInfo/${species}/chr/" . $transcriptInfo{$hitName}{$species}{'chr'};
      die "Error[extractFeatures.pl]! cannot find the $fileChr for the hitName $hitName, maybe different format $!\n" unless (-f "$fileChr");
      my $string;
      @{$transcriptInfo{$hitName}{$species}{'exons'}} = sort {$a->{"start"} <=> $b->{"start"} } @{$transcriptInfo{$hitName}{$species}{'exons'}};

      #concatenate fasta
      if ($concatenate eq "yes"){
	foreach my $exon (@{$transcriptInfo{$hitName}{$species}{'exons'}}){
	  my $command = "$chr_subseq $fileChr " . $exon->{"start"} . ' ' . $exon->{"end"};
	  $string .= `$command`; die ("Error[extractFeatures.pl]! chr_subseq didn't work for the command $command:$!\n") if ($?);
	  chomp $string;
	}
	print OUT ">$species\n";
	if ($transcriptInfo{$hitName}{$species}{'strand'} eq '-'){
	  $string = getReverseComplement($string);
	}
	$string = string2FASTA($string);
	chomp $string;
	print OUT "$string\n";
      }

      #do not concatenate fasta
      if($concatenate eq "no"){
	foreach my $exon (@{$transcriptInfo{$hitName}{$species}{'exons'}}){
	  my $command = "$chr_subseq $fileChr " . $exon->{"start"} . ' ' . $exon->{"end"};
	  $string = `$command`; die ("Error[extractFeatures.pl]! chr_subseq didn't work for the command $command:$!\n") if ($?);
	  chomp $string;
	  print OUT ">${species}_exon" . $exon->{"start"} . '_' . $exon->{"end"} . "\n";
	  if ($transcriptInfo{$hitName}{$species}{'strand'} eq '-'){
	    $string = getReverseComplement($string);
	  }
	  $string = string2FASTA($string);
	  chomp $string;
	  print OUT "$string\n";
	}
      }
    }
  close OUT;
  }
}






