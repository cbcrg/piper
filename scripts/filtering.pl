#!/usr/bin/env perl
use strict;
use warnings;
#Author: Giovanni Bussotti

#HELP
my $help  = 0;
my $cmd;
foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '-help') or ($ARGV[$field] eq '--help'));
}
if ($help > 0){
 help_message();
}


#TAKE OPTIONS
acceptedVariableSpace();
my ($clusterConfigLine , @clusterHeaderLines , @allSplitNames);
my ($querySplits4cluster4filtering , $uniprot , $overlapStrand , $geneidParameterFile , $rfam , $blastName , $orf_score_threshold , $geneid , $cluster , $codingPotential , $annotationFile , $overlapDistance , $nr , $pfam , $blastx , $rpsblast , $experimentName , $pipelineDirName) = options();
my $queryName          = "${pipelineDirName}/experiments/${experimentName}/pre_processing/inputInfo/___query___.fa";
my $gtfTxName          = "${pipelineDirName}/experiments/${experimentName}/pre_processing/inputInfo/___tx___.gtf";
my $gtfExName          = "${pipelineDirName}/experiments/${experimentName}/pre_processing/inputInfo/___ex___.gtf";
my $ncbiBlastXConfig   = "${pipelineDirName}/experiments/${experimentName}/CONFIG/ncbiBlastXConfig";
my $ncbiRPSblastConfig = "${pipelineDirName}/experiments/${experimentName}/CONFIG/ncbiRPSblastConfig";
my $geneidConfig       = "${pipelineDirName}/experiments/${experimentName}/CONFIG/geneidConfig";
my $blastConfig        = "${pipelineDirName}/experiments/${experimentName}/CONFIG/blast4rfamFiltering";
my $clusterConfigName  = "${pipelineDirName}/experiments/${experimentName}/CONFIG/clusterConfig";

exit if (-f "${pipelineDirName}/experiments/${experimentName}/pre_processing/inputInfo/filterRun");

open (F , "<$queryName") or die "Error[filtering.pl]! Cannot read the query file $queryName $!\n";
my %fasta = loadMultifastaIntoHash(<F>);
close F;
my %infoExGtf  = readGTF($gtfExName);
my %infoTxGtf  = readGTF($gtfTxName);
my ($ncbiBlastXConfigLine , $ncbiRPSblastConfigLine) = readNcbiConfig($ncbiBlastXConfig , $ncbiRPSblastConfig);



#OVERLAP WITH ANNOTATIONS
if (($annotationFile ne "none") && (scalar keys %fasta != 0)){
  print STDERR "#REMOVING TRANSCRIPTS OVERLAPPING WITH PROTEIN ANNOTATIONS..\n";
  my $s;
  if ($overlapStrand eq 'yes'){$s = 1;}
  else {$s = 0;}
  my  $gtfTxNameExt;
  if ($overlapDistance > 0){
    $gtfTxNameExt = extendRealGtf2File();
  }
  else {$gtfTxNameExt = $gtfTxName;}
  $cmd = "${pipelineDirName}/scripts/overlap $gtfTxNameExt $annotationFile -v -st $s 2> /dev/null";
  my @overlapOut = `$cmd`;  die "Error[filtering.pl]! Error with command:\n$cmd\n$!\n" if ($?);
  (system "rm $gtfTxNameExt") == 0 or die "Error[filtering.pl]! cannot remove $gtfTxNameExt\n$!\n" if ($overlapDistance > 0);
  foreach my $line (@overlapOut){
    my $transcriptName;
    if ($line=~/transcript_id \"([^\"]+)\"/){
      $transcriptName = $1;
    }
    else{die "Error[filtering.pl]! Cannot find the transcript name in overlap out line\n$line\n";}
    if ($line=~/ov_feat2: 1/){
      print "#$transcriptName\n";
      delete($fasta{$transcriptName});
      delete($infoExGtf{$transcriptName});
      delete($infoTxGtf{$transcriptName});
    }
    elsif ($line=~/ov_feat2: 0/){}
    else {die "Error[filtering.pl]! Cannot find the ov_feat2: flag in overlap out line\n$line\n";}
  }
  printRemaining();
}



#GENEID
if (($codingPotential eq 'on') && (scalar keys %fasta != 0)){
  print STDERR "#REMOVING TRANSCRIPTS HAVING GENEID CODING POTENTIAL..\n";
  my $geneIdOptions;
  open (GIC,"$geneidConfig") or die "Error[filtering.pl]! Cannot read the geneidConfig file\n$!\n";
  foreach my $line (<GIC>){
    chomp $line;
    $geneIdOptions .= $line;
  }
  close GIC;
  die "Error[filtering.pl]! Please put back the \"-soW\" in the geneidConfig\n" unless ($geneIdOptions=~/-soW/);
  my @geneIdOut = `geneid -P $geneidParameterFile $geneIdOptions $queryName`;  die ("Error[filtering.pl]! geneid didn't work for the command \n\'geneid $geneIdOptions $queryName\'\n$!\n") if ($?);

  my %h_dat = parsingGeneIdOut (@geneIdOut);
  foreach my $transcriptName (keys %h_dat) {
    $h_dat{$transcriptName}->{'orf_score'} = 0 unless (defined $h_dat{$transcriptName}->{'orf_score'});
    if ($h_dat{$transcriptName}->{'orf_score'} > $orf_score_threshold){
      print "#$transcriptName\n";
      delete($fasta{$transcriptName});
      delete($infoExGtf{$transcriptName});
      delete($infoTxGtf{$transcriptName});
    }
  }
  printRemaining();
}


#RFAM
if (($rfam ne "none") && (scalar keys %fasta != 0)){
  print STDERR "#REMOVING TRANSCRIPTS HAVING MATCHING AN RFAM TRANSCRIPT..\n";
  my $configLine = readBlastConfig ($blastConfig);
  $cmd = "$blastName $rfam  $queryName -mformat=2 $configLine";
  my @blastOUT = `$cmd`;  die ("Error[filtering.pl]! blast didn't work for the command \'$cmd\'\n$!\n") if ($?);
  my %alreadySeen;
  if (@blastOUT){
    foreach my $line (@blastOUT){
      chomp $line;
      if ($line=~/^(\S+)/){
	my $transcriptName = $1;
	print "#$transcriptName\n" if (! defined $alreadySeen{$transcriptName});
	delete($fasta{$transcriptName});
	delete($infoExGtf{$transcriptName});
	delete($infoTxGtf{$transcriptName});
	$alreadySeen{$transcriptName} = 1;
      }
    }
    printRemaining();
  }
}




#RPS-BLAST
if (($pfam ne "none") && (scalar keys %fasta != 0)){
  print STDERR "#REMOVING TRANSCRIPTS HAVING A TRANSLATIONAL PRODUCT MATCHING A PFAM MODEL..\n";
  $cmd = "$rpsblast -db $pfam -query $queryName -outfmt 6 $ncbiRPSblastConfigLine";
  my @blastOUT = `$cmd`;  die ("Error[filtering.pl]! RPS-blast didn't work for the command \'$cmd\'\n$!\n") if ($?);
  my %alreadySeen;
  if (@blastOUT){
    foreach my $line (@blastOUT){
      chomp $line;
      if ($line=~/^(\S+)/){
	my $transcriptName = $1;
	print "#$transcriptName\n" if (! defined $alreadySeen{$transcriptName});
	delete($fasta{$transcriptName});
	delete($infoExGtf{$transcriptName});
	delete($infoTxGtf{$transcriptName});
	$alreadySeen{$transcriptName} = 1;
      }
    }
    printRemaining();
  }
}


#BLASTX vs UNIPROT
if (($uniprot ne "none") && (scalar keys %fasta != 0)){
  print STDERR "#REMOVING TRANSCRIPTS HAVING A TRANSLATIONAL PRODUCT MATCHING A UNIPROT PROTEIN..\n";
  my @blastOUT;
  if ($cluster eq 'on'){
    @allSplitNames = ();
    @blastOUT = runBlastXonTheCluster("uniprot");
  }
  else{
    $cmd = "$blastx -db $uniprot -query $queryName -outfmt 6  $ncbiBlastXConfigLine";
    @blastOUT = `$cmd`;  die ("Error[filtering.pl]! blastx didn't work for the command \'$cmd\'\n$!\n") if ($?);
  }
  my %alreadySeen;
  if (@blastOUT){
    foreach my $line (@blastOUT){
      chomp $line;
      if ($line=~/^(\S+)/){
	my $transcriptName = $1;
	print "#$transcriptName\n" if (! defined $alreadySeen{$transcriptName});
	delete($fasta{$transcriptName});
	delete($infoExGtf{$transcriptName});
	delete($infoTxGtf{$transcriptName});
	$alreadySeen{$transcriptName} = 1;
      }
    }
    printRemaining();
  }
}

#BLASTX vs NR
if (($nr ne "none") && (scalar keys %fasta != 0)){
  print STDERR "#REMOVING TRANSCRIPTS HAVING A TRANSLATIONAL PRODUCT MATCHING A NR PROTEIN..\n";
  my @blastOUT;
  if ($cluster eq 'on'){
    @allSplitNames = ();
    @blastOUT = runBlastXonTheCluster("nr");
  }
  else{
    $cmd = "$blastx -db $nr -query $queryName -outfmt 6  $ncbiBlastXConfigLine";
    @blastOUT = `$cmd`;  die ("Error[filtering.pl]! blastx didn't work for the command \'$cmd\'\n$!\n") if ($?);
  }
  my %alreadySeen;
  if (@blastOUT){
    foreach my $line (@blastOUT){
      chomp $line;
      if ($line=~/^(\S+)/){
	my $transcriptName = $1;
	print "#$transcriptName\n" if (! defined $alreadySeen{$transcriptName});
	delete($fasta{$transcriptName});
	delete($infoExGtf{$transcriptName});
	delete($infoTxGtf{$transcriptName});
	$alreadySeen{$transcriptName} = 1;
      }
    }
    printRemaining();
  }
}





system "touch ${pipelineDirName}/experiments/${experimentName}/pre_processing/inputInfo/filterRun";






####FUNCTIONS
sub options {
  my ($querySplits4cluster4filtering , $uniprot , $overlapStrand , $geneidParameterFile , $rfam , $blastName , $orf_score_threshold , $geneid , $cluster , $codingPotential , $annotationFile , $overlapDistance , $nr , $pfam , $blastx , $rpsblast , $experimentName , $pipelineDirName);
  my $spyOrf_score_threshold           = 1;
  my $spyCluster                       = 1;
  my $spyAnnotationFile                = 1;
  my $spyOverlapDistance               = 1;
  my $spyNr                            = 1;
  my $spyUniprot                       = 1;
  my $spyPfam                          = 1;
  my $spyBlastX                        = 1;
  my $spyRpsBlast                      = 1;
  my $spyExperiment                    = 1;
  my $spyPipelineDir                   = 1;
  my $spyCodingPotential               = 1;
  my $spyGeneid                        = 1;
  my $spyBlast                         = 1;
  my $spyRfam                          = 1;
  my $spyGeneidParameterFile           = 1;
  my $spyOverlapStrand                 = 1;
  my $spyQuerySplits4cluster4filtering = 1;

  foreach my $field (0..$#ARGV){
    if ($ARGV[$field] eq '-querySplits4cluster4filtering'){
	$querySplits4cluster4filtering = $ARGV[1+$field];
	$spyQuerySplits4cluster4filtering = 2;
	next;
    }
    if ($spyQuerySplits4cluster4filtering == 2){
      $spyQuerySplits4cluster4filtering = 3;
      next;
    }
    if ($ARGV[$field] eq '-overlapStrand'){
	$overlapStrand = $ARGV[1+$field];
	$spyOverlapStrand = 2;
	next;
    }
    if ($spyOverlapStrand == 2){
      $spyOverlapStrand = 3;
      next;
    }
    if ($ARGV[$field] eq '-rfam'){
      $rfam = $ARGV[1+$field];
      $spyRfam = 2;
      next;
    }
    if ($spyRfam == 2){
      $spyRfam = 3;
      next;
    }
    if ($ARGV[$field] eq '-blast'){
      $blastName = $ARGV[1+$field];
      $spyBlast = 2;
      next;
    }
    if ($spyBlast == 2){
      $spyBlast = 3;
      next;
    }
    if ($ARGV[$field] eq '-orf_score_threshold'){
	$orf_score_threshold = $ARGV[1+$field];
	$spyOrf_score_threshold = 2;
	next;
    }
    if ($spyOrf_score_threshold == 2){
	$spyOrf_score_threshold = 3;
	next;
    }
    if ($ARGV[$field] eq '-geneid'){
	$geneid = $ARGV[1+$field];
	$spyGeneid = 2;
	next;
    }
    if ($spyGeneid == 2){
	$spyGeneid = 3;
	next;
    }
    if ($ARGV[$field] eq '-codingPotential_check'){
	$codingPotential = $ARGV[1+$field];
	$spyCodingPotential = 2;
	next;
    }
    if ($spyCodingPotential == 2){
	$spyCodingPotential = 3;
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
    if ($ARGV[$field] eq '-rpsblast'){
	$rpsblast = $ARGV[1+$field];
	$spyRpsBlast = 2;
	next;
    }
    if ($spyRpsBlast == 2){
      $spyRpsBlast = 3;
      next;
    }
    if ($ARGV[$field] eq '-blastx'){
	$blastx = $ARGV[1+$field];
	$spyBlastX = 2;
	next;
    }
    if ($spyBlastX == 2){
      $spyBlastX = 3;
      next;
    }
    if ($ARGV[$field] eq '-pfam'){
	$pfam = $ARGV[1+$field];
	$spyPfam = 2;
	next;
    }
    if ($spyPfam == 2){
      $spyPfam = 3;
      next;
    }
    if ($ARGV[$field] eq '-nr'){
	$nr = $ARGV[1+$field];
	$spyNr = 2;
	next;
    }
    if ($spyNr == 2){
      $spyNr = 3;
      next;
    }
    if ($ARGV[$field] eq '-uniprot'){
	$uniprot = $ARGV[1+$field];
	$spyUniprot = 2;
	next;
    }
    if ($spyUniprot == 2){
      $spyUniprot = 3;
      next;
    }
    if ($ARGV[$field] eq '-overlapDistance'){
	$overlapDistance = $ARGV[1+$field];
	$spyOverlapDistance = 2;
	next;
    }
    if ($spyOverlapDistance == 2){
      $spyOverlapDistance = 3;
      next;
    }
    if ($ARGV[$field] eq '-annotation'){
	$annotationFile = $ARGV[1+$field];
	$spyAnnotationFile = 2;
	next;
    }
    if ($spyAnnotationFile == 2){
      $spyAnnotationFile = 3;
      next;
    }
    if ($ARGV[$field] eq '-cluster'){
	$cluster = $ARGV[1+$field];
	$spyCluster = 2;
	next;
    }
    if ($spyCluster == 2){
	$spyCluster = 3;
	next;
    }
    if ($ARGV[$field] eq '-geneid_parameter'){
	$geneidParameterFile = $ARGV[1+$field];
	$spyGeneidParameterFile = 2;
	next;
    }
    if ($spyGeneidParameterFile == 2){
	$spyGeneidParameterFile = 3;
	next;
    }
  }
  if ((defined $overlapStrand) && ($overlapStrand ne 'yes') && ($overlapStrand ne 'no')){
    die "Error[startPipeline.pl]! The -overlapStrand parameter accepts either \"yes\" or \"no\"\n";
  }
  if ((defined $codingPotential) && ($codingPotential ne 'on') && ($codingPotential ne 'off')){
    die "Error[startPipeline.pl]! The -codingPotential_check parameter accepts either \"on\" or \"off\"\n";
  }
  $cluster             = 'off'                 if (! defined $cluster);
  $blastx              = 'blastx'              if (! defined $blastx);
  $rpsblast            = 'rpsblast'            if (! defined $rpsblast);
  $annotationFile      = "none"                if (! defined $annotationFile);
  $overlapStrand       = "yes"                 if (! defined $overlapStrand);
  $overlapDistance     = 0                     if (! defined $overlapDistance);
  $nr                  = "none"                if (! defined $nr);
  $uniprot             = "none"                if (! defined $uniprot);
  $pfam                = "none"                if (! defined $pfam);
  $geneid              = "geneid"              if (! defined $geneid);
  $codingPotential     = "off"                 if (! defined $codingPotential);
  $orf_score_threshold = 20                    if (! defined $orf_score_threshold);
  $rfam                = 'none'                if (! defined $rfam);
  $blastName           = 'wu-blastn'           if (! defined $blastName);
  $geneidParameterFile = 'none'                if (! defined $geneidParameterFile);
  $querySplits4cluster4filtering = 25          if (! defined $querySplits4cluster4filtering);
  die "Error[filtering.pl]! You must provide the -experiment parameter.\n"   if ($spyExperiment  != 3);
  die "Error[filtering.pl]! You must provide the -pipeline_dir parameter.\n" if ($spyPipelineDir != 3);
  die "Error[filtering.pl]! If you want to use geneid you must specify with -geneid_parameter the parameter file estimated on your reference species. You can find a list of precomputed parameter files here: http://genome.crg.es/software/geneid/index.html#parameters either you can estimate a new one: http://genome.crg.es/software/geneid/training.html\n" if (($geneidParameterFile eq 'none') && ($codingPotential eq 'on'));

  return ($querySplits4cluster4filtering , $uniprot , $overlapStrand , $geneidParameterFile , $rfam , $blastName , $orf_score_threshold , $geneid , $cluster , $codingPotential , $annotationFile , $overlapDistance , $nr , $pfam , $blastx , $rpsblast , $experimentName , $pipelineDirName);
}

sub acceptedVariableSpace {
  my %space = ('-querySplits4cluster4filtering' => 1 , '-uniprot' => 1 , '-overlapStrand' => 1  , '-orf_score_threshold' => 1 , '-rfam' => 1  , '-blast' => 1  , '-codingPotential_check' => 1  , '-geneid' => 1  , '-experiment' => 1  , '-pipeline_dir' => 1  , '-rpsblast' => 1  , '-blastx' => 1  , '-overlapDistance' => 1  , '-pfam' => 1  , '-nr' => 1  , '-geneid_parameter' => 1  , '-annotation' => 1  , '-cluster' => 1 );
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[filtering.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}

sub help_message {
my $helpMessage = "\nNAME
filtering.pl - It accepts the ___ex___.gtf    ___query___.fa  ___tx___.gtf files created  by createInput.pl in the \"pre_processing/inputInfo/\" input folder. It runs several filters to get a putative ncRNA set

SYNOPSIS
filtering.pl -experiment -pipeline_dir [-orf_score_threshold -rfam -blast -codingPotential_check -geneid  -rpsblast -blastx -overlapDistance -pfam -nr -geneid_parameter -annotation -cluster]

DESCRIPTION
   * It runs the following filtering steps:
   * transcripts VS annotations (overlap)
   * transcripts VS coding potential check (geneid)
   * transcript  VS rfam (blast)
   * transcript  VS pfam (rpsblast)
   * transcripts VS nr (blastx)

OPTIONS
   ***ANNOTATION FILTER
      * To activate this filter you have to provide a gtf annotation file via -annotation parameter
      * By default the overlap distance to the annotations is 0. You can inscrease via -overlapDistance
      * By default the overlap is given to features laying on the same strand (-overlapStrand \"yes\"). You can remove a feature if it is overlapping on the other strand by doing -overlapStrand \"no\"

   ***CODING POTENTIAL FILTER
      * To activate this filter you have to set \"-codingPotential_check on\" and you have to specify the geneid parameter file of your reference species via -geneid_parameter
      * By default geneid program name is \"geneid\". If this is not the case you must provide the field -geneid with the -geneid path.

   ***RFAM FILTER
      * To activate this filter you have to specify via -rfam the xdformat formatted Rfam database
      * By default blast program name is \"wu-blast\". If this is not the case you must provide the field -blast with the blast path.

   ***PFAM FILTER
      * To activate this filter you have to specify via -pfam the Pfam database formatted for rpsblast
      * You should be able to find Pfam already formatted here ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/
      * By default rpsblast program name is \"rpsblast\". If this is not the case you must provide the field -rpsblast with the rpsblast path.

   ***BLASTX FILTER
      * To activate this filter you have to specify via -nr the nr the formatdb formatted nr database, and/or -uniprot the uniprot formatdb formatted database
      * By default blastx program name is \"blastx\". If this is not the case you must provide the field -blastx with the blastx path.

TROUBLESHOOTING
By default no filter is applyied! The user must select the filters he wants to apply
Once the script finished it will create a file called \"/pre_processing/inputInfo/filterRun\". The user must remove this file in the case he wants to rerun filtering.pl
this filtering routine accepts ncbi-blastx+ 

";
print "$helpMessage\n\n\n";
exit;
}

sub printRemaining {
  open (F , ">$queryName") or die "Error[filtering.pl]! Cannot re-write the $queryName file\n";
  open (E , ">$gtfExName") or die "Error[filtering.pl]! Cannot re-write the $gtfExName file\n";
  open (T , ">$gtfTxName") or die "Error[filtering.pl]! Cannot re-write the $gtfTxName file\n";
  if (scalar keys %infoTxGtf == 0 ){print STDERR "#All the queries have been removed! Exiting\n"; exit;}
  foreach my $remainingTranscript (keys %fasta){
    print F ">$remainingTranscript\n$fasta{$remainingTranscript}\n";
    print E $infoExGtf{$remainingTranscript} . "\n";
    print T $infoTxGtf{$remainingTranscript} . "\n";
  }
  close F;
  close E;
  close T;
}


sub loadMultifastaIntoHash{
  my (@multifasta) = @_;
  my (%allTheSequences , $sequence , $header);
  my $spy = 1;
  foreach my $line (@multifasta){
    chomp $line;
    next if ($line=~/^\s*$/);
    if ($line=~/>(\S+)/){
      if ($spy == 2){
        $spy = 1;
        $allTheSequences{$header} = $sequence;
        $header   = '';
        $sequence = '';
      }
      $header = $1;
      next;
    }
    if (($line!~/>/) and ($line=~/\w+/)){
      $sequence .= $line;
      $spy = 2;
      next;
    }
  }
  $allTheSequences{$header} = $sequence;
  return %allTheSequences;
}
sub trim {
  my ($string) = @_;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
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

sub readGTF {
  my ($gtfFile) = @_;
  my %hash;
  open (G , "<$gtfFile") or die "Error[filtering.pl]! Cannot read $gtfFile $!\n"; 
  foreach my $line (<G>){
    chomp $line;
    my $transcriptName;
    if ($line=~/transcript_id \"([^\"]+)\"/){
       $transcriptName = $1;
       $hash{$transcriptName} = $line;
    }
    else{
      die "Error[filtering.pl]! Cannot find the transcript name in $gtfFile\nline\n$line\n$!\n";
    }
  }
  close G;
  return %hash;
}

sub extendRealGtf2File {
  open (EXT,">${gtfTxName}_ext") or die "Error[filtering.pl]! Cannot create the ${gtfTxName}_ext file $!\n";
  open (G2,"<$gtfTxName") or die "Error[filtering.pl]! cannot open the -gtf2 file $!\n";
  foreach my $line (<G2>){
    if ($line=~/^(\S+)(\s+\S+\s+\S+\s+)(\S+)\s+(\S+)(\s+.+)/){
      my $chr           = $1;
      my $block1        = $2;
      my $extendedStart = min2 ($3 , $4) - $overlapDistance;
      my $extendedEnd   = max2 ($4 , $3) + $overlapDistance;
      my $block2        = $5;
      $extendedStart    = 1 if ($extendedStart < 0);
      ##$extendedEnd      = $infoAssembly{$chr} if ($extendedEnd > $infoAssembly{$chr}) ;
      print EXT "$chr$block1$extendedStart\t$extendedEnd$block2\n";
    }
  }
  close G2;
  close EXT;
  return "${gtfTxName}_ext";
}

sub parsingGeneIdOut {
  my @geneIdOut = @_;
  my ($id , %h_dat);
  foreach my $line (@geneIdOut) {
    chomp ($line);
    next if ($line =~ /^##/);
    if ($line =~ /^# Sequence (\S+) - Length = (\d+) bps$/) {
      $id              = $1;
      my $length_tr    = $2;
      $h_dat{$id}->{'tr_length'} = $length_tr;
    }
    elsif ($line =~ /^\s+Single\s+.*/) {
      my @tab= split (/\s+/, $line);
      my $score           = $tab[10];
      my $orf_length      = $tab[12];
      $h_dat{$id}->{'orf_length'}  = max2($h_dat{$id}->{'orf_length'}, $orf_length);
      $h_dat{$id}->{'orf_score'}   = max2($h_dat{$id}->{'orf_score'}, $score);
    }
    elsif ($line =~ /^# Singles(.*)/){
    }
    else {
      print "\n#Error: could not parse line:\n".$line."\n";
      exit 1;
    }
  }
  return %h_dat;
}

sub readBlastConfig {
  my ($configName) = @_;
  my $configLine;
  open (C,"<$configName") or die "Error[filtering.pl]! Cannot read the blast configuration file $!\n";
  foreach my $line (<C>){
    chomp $line;
    next if ($line=~/^\s*$/);
    $configLine .= $line . " ";
  }
  if ($configLine =~/mformat/){ die "Error[filtering.pl]! Please erease from the blastConfig the mformat option \n";}
  close C;
  return $configLine;
}

sub fileNameGenerator{
  my ($nameRoot) = @_;
  my $tmp_name_counter = 0;
  my $tmp_name;
  while (!$tmp_name || -f $tmp_name) {
    $tmp_name_counter++;
    $tmp_name = "${nameRoot}_$$".".${tmp_name_counter}.fa";
  }
  return $tmp_name;
}

sub readNcbiConfig {
  my ($configName1 , $configName2) = @_;
  my ($ncbiBlastXConfigLine , $ncbiRPSblastConfigLine);
  open (C,"<$configName1") or die "Error[filtering.pl]! Cannot read the ncbi blastX configuration file $!\n";
  foreach my $line (<C>){
    chomp $line;
    next if ($line=~/^\s*$/);
    $ncbiBlastXConfigLine .= $line . " ";
  }
  close C;
  if ( $ncbiBlastXConfigLine =~/-outfmt/){ die "Error[filtering.pl]! Please erease from the ncbiBlastXConfig the -outfmt option \n";}

  open (C,"<$configName2") or die "Error[filtering.pl]! Cannot read the ncbi rpsBlast configuration file $!\n";
  foreach my $line (<C>){
    chomp $line;
    next if ($line=~/^\s*$/);
    $ncbiRPSblastConfigLine .= $line . " ";
  }
  close C;
  if ( $ncbiRPSblastConfigLine =~/-outfmt/){ die "Error[filtering.pl]! Please erease from the ncbiRPSblastConfig the -outfmt option \n";}

  return ($ncbiBlastXConfigLine , $ncbiRPSblastConfigLine);
}

sub waiting {
  my @allSplitNames= @_;
  my $total     = scalar @allSplitNames ;
  my $time = 15;
  my $done = 0;
  while ($done < $total ){
    print STDERR "#Waiting the job to exit the cluster... " . $done . " already finished\n";
    system "sleep $time";
    #$time += 20;
    $done = `ls ${pipelineDirName}/experiments/${experimentName}/pre_processing/inputInfo/experiments/inputInfo/CLUSTER_FILES/*___blastx_done___ 2> /dev/null | wc -l`;
    chomp $done;
  }
  print STDERR "#All jobs exited the cluster\n";
}
sub readClusterConfig {
    open (CC,"<$clusterConfigName") or die "Error! Cannot read $clusterConfigName \n$!\n";
    my $spyHeader = 0;
    foreach my $line (<CC>){
	chomp $line;
	next if ($line=~/^\s*$/);
	$clusterConfigLine = $line if (($line=~/##SCRIPT##/) and ($spyHeader == 0));

	if ($line=~/^#HEADER/){$spyHeader = 1;next;}
	push (@clusterHeaderLines , $line) if ($spyHeader == 1);
    }
    close CC;
}
sub clearDone {
    my $check_blast_done_ =  `ls ${pipelineDirName}/experiments/${experimentName}/pre_processing/inputInfo/experiments/inputInfo/CLUSTER_FILES/*___blastx_done___ 2> /dev/null | wc -l`;
    chomp $check_blast_done_;
    if ($check_blast_done_ > 0){
      my $rmCmd = "rm ${pipelineDirName}/experiments/${experimentName}/pre_processing/inputInfo/experiments/inputInfo/CLUSTER_FILES/*___blastx_done___";
      (system "$rmCmd")==0 or die "Error [filtering.pl]! with command $rmCmd \n$!\n";
    }
}

sub splitTheQueries {
  my ($targetDBname) = @_;
  #split the query file
  my $splits = $querySplits4cluster4filtering;
  my $tmp_split_dir       = "${pipelineDirName}/experiments/${experimentName}/pre_processing/inputInfo/experiments/inputInfo/CLUSTER_FILES";
  my $queriesForEachSplit = (scalar keys %fasta) / ($splits);
  $queriesForEachSplit    = sprintf("%.0f", $queriesForEachSplit);
  $queriesForEachSplit = 0 if ($queriesForEachSplit == 0);

  #fileNameGenerator
  my $c = 0;
  my $tmp_split_name = fileNameGenerator("${tmp_split_dir}/split4nrSearch_$targetDBname");
  push (@allSplitNames , $tmp_split_name);
  open (TF,">>$tmp_split_name") or die "Error [filtering.pl]! Cannot create $tmp_split_name in the runBlastXonTheCluster function\n$!\n";
  foreach my $query (keys %fasta){
    $c++;
    if ($c <= $queriesForEachSplit){
      print TF ">$query\n$fasta{$query}\n";
    }
    else {
      close TF;
      $c = 1;
      $tmp_split_name = fileNameGenerator("${tmp_split_dir}/split4nrSearch_$targetDBname");
      push (@allSplitNames , $tmp_split_name);
      open (TF,">>$tmp_split_name") or die "Error [filtering.pl]! Cannot create $tmp_split_name in the runBlastXonTheCluster function\n$!\n";
      print TF ">$query\n$fasta{$query}\n";
    }
  }
  close TF;
}


sub runBlastXonTheCluster {
  my ($targetDBname) = @_;
  print STDERR "#run blastX on the cluster..\n";

  splitTheQueries("$targetDBname");

  #create cluster scripts for each split
  readClusterConfig();
  foreach my $split (@allSplitNames){
    my $scriptName = "${split}.sh";
    open (S,">$scriptName") or die "Error [filtering.pl]! Cannot create $scriptName in the runBlastXonTheCluster function\n$!\n";
    foreach my $h (@clusterHeaderLines){
      print S "$h\n";
    }
    my $blastx_cmd ;
    if ($targetDBname eq 'nr'){
      $blastx_cmd = "$blastx -db $nr -query $split -outfmt 6 $ncbiBlastXConfigLine >> ${split}.o 2>> ${split}.e";
    }
    if ($targetDBname eq 'uniprot'){
      $blastx_cmd = "$blastx -db $uniprot -query $split -outfmt 6 $ncbiBlastXConfigLine >> ${split}.o 2>> ${split}.e";
    }
    print S "$blastx_cmd" . "\n\n";
    print S "touch ${split}.___blastx_done___" . "\n";
    close S;

    #do the cmd
    my $cmd = $clusterConfigLine;
    $cmd =~s/##SCRIPT##/$scriptName/;
    $cmd =~s/-o \S+/-o \/dev\/null/;
    $cmd =~s/-e \S+/-e \/dev\/null/;
    (system "$cmd") == 0 or die "Error[filtering.pl]! cannot run on the cluster:\n$cmd\n$!\n";
  }

  waiting(@allSplitNames);
  clearDone();
  clusterErrorCheck();

  #read the outputs
  my @toRemove;
  foreach my $split (@allSplitNames){
    my $outName = "${split}.o";
    open (O,"<$outName") or die "Error [filtering.pl]! Cannot read $outName blastx output \n$!\n";
    foreach my $line (<O>){
      if ($line=~/^(\S+)/){
	push (@toRemove , $1);
      }
    }
    close O;
  }
  return @toRemove;
}

sub clusterErrorCheck {
  my $tmp_split_dir       = "${pipelineDirName}/experiments/${experimentName}/pre_processing/inputInfo/experiments/inputInfo/CLUSTER_FILES";
  my @filesHavingError;
  my @e1 = `grep -l alloc ${tmp_split_dir}/*.e`;
  my @e2 = `grep -l ERROR ${tmp_split_dir}/*.e`;
  my @e3 = `grep -l error ${tmp_split_dir}/*.e`;
  my @e4 = `grep -l Error ${tmp_split_dir}/*.e`;
#  my @e5 = `grep -l FATAL ${clusterFileDir}/*`;
  my @e6 = `grep -l Killed ${tmp_split_dir}/*.e`;
  my @e7 = `grep -l memory ${tmp_split_dir}/*.e`;


  push (@filesHavingError , @e1 ) if (@e1);
  push (@filesHavingError , @e2 ) if (@e2);
  push (@filesHavingError , @e3 ) if (@e3);
  push (@filesHavingError , @e4 ) if (@e4);
 # push (@filesHavingError , @e5 ) if (@e5);
  push (@filesHavingError , @e6 ) if (@e6);
  push (@filesHavingError , @e7 ) if (@e7);


  if (@filesHavingError){
    print "Error[filtering.pl]! The following files contain the word \"error\" or \"malloc\"  or \"Killed\" or \"memory\" while running on the cluster:\n";
    foreach my $f (@filesHavingError){
      print $f ;
    }
    die "Please rerun manually the individual shell scripts before listed\n";
  }
}
