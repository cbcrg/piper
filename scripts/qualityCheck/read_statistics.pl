#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use File::Basename;
#Author: Giovanni Bussotti

#HELP
my $help  = 0;
foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '-help') or ($ARGV[$field] eq '--help'));
}
if ($help > 0){
 help_message();
}

#GLOBAL VARIABLES
my ($exclude , $type , $repeatThreshold , $splitMode , $randomizations , $shuffleBed , $shuffleBed_gap , $millionMappedReads , $tmp_split_name , $exonGtfSplits , $cluster , $species , $stranded , $bamtools , $bamdir ,   $experimentName , $pipelineDirName) = options();
acceptedVariableSpace();
$randomizations++ if ($randomizations > 0);
my $bamdir_name = basename($bamdir);
my $outDir      = "${pipelineDirName}/experiments/${experimentName}/qualityCheck/reads_statistics/${bamdir_name}_VS_${species}_${type}_rep$repeatThreshold/";
my ($tx_gtf , $gtf_map , %tx_info , %ex_info , @allExons , @allSplitNames , $clusterConfigLine ,@clusterHeaderLines );
$gtf_map = set_gtf_map();
my $gtfMapSaveName = $gtf_map;

my $count                = 0;
(system "rm -rf $outDir") == 0 or die "Error[reads_statistics.pl]! cannot remove the directory $outDir  \n$!\n" if ((-d $outDir)&&($splitMode eq 'off')) ;
(system "mkdir -p $outDir") == 0 or die "Error[reads_statistics.pl]! cannot create the $outDir directory \n$!\n";
printLog($outDir);
#BAM INDICES
checkIndices();
#COUNT ALL THE MAPPING READS
$millionMappedReads = countAllMappingReads() if (! defined $millionMappedReads);
#TAKE ASSEMBLY INFO
my ($shuffleBed_chrSize , %gtf_map_info_ex , %gtf_map_info_tx , $gtf_map_block , $numberTranscripts , $syntaxCheck , @deltas4block , @blocks , @blockStarts);
#CHECK SYNTAX BAM FILE
$syntaxCheck = syntax_check();

if ($randomizations > 0) {
    my $assemblyDir       = "${pipelineDirName}/allGenomeInfo/${species}/chr/";
    my %infoAssembly      = takeAssemblyInfo($assemblyDir);
    $shuffleBed_chrSize   = makeShufflebedChrSize(%infoAssembly);
    %gtf_map_info_ex      = readingGTFexons($gtf_map);
    %gtf_map_info_tx      = gtfExons2Transcripts(%gtf_map_info_ex);
    takeBlocks($gtf_map);
    $gtf_map_block        = doBlocks();
}

while ( $randomizations >= 0 ){
    $randomizations--;
    my $ran = fileNameGenerator ('ran_');
    initialize();

    $gtf_map = $gtfMapSaveName;
    if (($splitMode eq 'off')&&($randomizations > 0)) {createRandomWithShuffleBed($ran);  $gtf_map = $ran; print STDERR "#using a gtf randomization..\n";}
    if (($splitMode eq 'off')&&($randomizations == 0)) {print STDERR "#using $gtf_map ..\n";}
    $gtf_map    = $tmp_split_name if ((defined $tmp_split_name) && ($splitMode eq 'on'));

    if ($cluster eq 'on'){
	my $splitDir = "${pipelineDirName}/experiments/${experimentName}/CLUSTER_FILES/splittedExons/";
	splitExons($splitDir);
	runOnTheCluster($splitDir);
	system "rm $gtfMapSaveName"     if ((defined $exclude) && ($gtfMapSaveName =~/^gtf_map_exclude_tmp.*/) && (-e $gtfMapSaveName));
	waiting($splitDir);
	concatenate();
	takeAverage();
	system "rm $ran" if (-f $ran);
	if ($randomizations > 0){
	  system "rm ${outDir}/exons.rpkm ; rm ${outDir}/transcripts.rpkm";
	  next;
	}
	else {
	    last;
	}
    }
    #COUNT READS PER EXONS
    countReadsPerExon();
    #OUTPUT
    my $tag = basename($gtf_map);
    printOut ($tag);
    if ($splitMode eq 'off'){
	renaming()   ;
	takeAverage();
        if ($randomizations > 0){system "rm ${outDir}/exons.rpkm ; rm ${outDir}/transcripts.rpkm";}
	system "rm $ran" if (-f $ran);
    }
    last if ($randomizations == 0);
}
system "rm $shuffleBed_chrSize" if (defined $shuffleBed_chrSize);
system "rm $gtf_map_block"      if (defined $gtf_map_block);
system "rm $gtfMapSaveName"     if ((defined $exclude) && ($gtfMapSaveName =~/^gtf_map_exclude_tmp.*/) && (-e $gtfMapSaveName));












sub help_message {
my $helpMessage = "\nNAME
read_statistics.pl

SYNOPSIS
perl read_statistics.pl -species -bamdir -experiment -pipeline_dir -type -repeatThreshold [-exclude -randomizations -shuffleBed -shuffleBed_gap -tmp_split_name -exonGtfSplits -cluster -bamtools -stranded  -millionMappedReads -help]

DESCRIPTION

- This script generate exon and transcript rpkm distributions given an input directory containing mapped reads (.bam). You can also choose to randomize the exons with shuffleBed to have the statistical significance of the overlap. This script will return also an average.tx and average.ex files showing the average rpkm values observed in each of the randomization. The real average values will be the bottom line of each file. The randomization strategy is prety complex and includes:
1- considering the input transcripts, any time the exon of one transcript overlap with the exon of another, the two transcript will be joined in the same \"block\". So each block is a collection of overlapping transcripts (that overlap at exonic level).
2- randomly map with shuffleBed all the blocks
3- split the projectd block back to exon annotations
This randomization strategy has the merit to mantain exactly the same nucleotide coverage as the original. Moreover the same exonic relative position is mantained. This means that the mapping of the block make it possible to mantain exactly the same amount of redundacy. In other words, the randomized exons at the end will be overlapping one another exactly as the original exons.
Moreover this script considers the homologs returned by pipeR, where each different hitName is a different transcript homolog. This important because the script is going to work both when pipeR returns for each query a single best hit (a single homolog), and when pipeR it is run exahustively, i.e. returning a one2many mapping. In the case pipeR generated multiple homologs for one query, each homolog will be considered as a different transcript, with a different id. This script will tell you how much is the average RPKM support of all the annotations, with respect to random projections of such annotations. The user can choose the repeatThreshold and the type (bh or all) of pipeR mapping.

OPTIONS
  *-species = <name>
  one of the species name found in the pipelineDir/allGenomeInfo   (human, for instance)

  *-bambir = <dirName>
  directory containing all the bam files of reads mapped on the same assembly as the one specified by species

  *-experiment = <name>
  experiment name (the species.repXX.ex.gtf will be taken by exprerimentName/results)

  *-pipeline_dir = <dir>
  pipeline position

  *-repeatThreshold = <int>
  choose as gtf_map the repeat threshold from the results folder. There are always 100 80 and 80. If the user run pipeR with another repeatThreshold, he can use it

  *-type = <bh|all>
  choose as gtf_map the best hit (bh) or the all file.

  *-tmp_split_name = <file>
  OPTION: consider this exon gtf file insted the one found at exprerimentName/results

  *-exonGtfSplits = <integer>
  OPTION: Number (approximate) of the splits you wanna create to run on the cluster. The number is approximate because the split criterion depends also on the transcript models, i.e. the splits will always include all the exons belonging to a transcript. Default [80]

  *-cluster = <on|off>
  OPTION: activate the cluster mode. Default [off]

  *-bamtools = <name>
  OPTION: specify the bamtools name on your computer. Default [bamtools]

  *-stranded = <0|1>
  OPTION: set it to 1 if you have stranded reads. The default is [0], this means that it will count the reads mapping on both strands on that region

  *-shuffleBed = <name>
  OPTION: shuffleBed name to run the randomizations. Default [shuffleBed]

  *-shuffleBed_gap = <name>
  OPTION: define the gap file for shuffleBed if available (eterochromatine, centromere, telomere, assembly gaps...available from UCSC tables for some genome assemblies).

  *-randomizations = <integer>
  OPTION: set the number of randomizations to run. Default [0].

  *-millionMappedReads = <integer>
  OPTION: to skipp the million of mappe reads calculation (to have the RPKM normalization) if this number is already known

  *-exclude = <string>
  OPTION: use this parameter to get rid of gtf_map lines including the specified string. Such string won't be consideted in the analysis. Use this option when your bam file do not contain a region that your gtf_map do.
  TIP: to check the chromosomes and their sizes, query the bam files with: samtools idxstats file.bam

";
print "$helpMessage\n\n\n";
exit;
}

sub initialize {
  %tx_info            = ();
  %ex_info            = ();
  @allExons           = ();
  @allSplitNames      = ();
  @clusterHeaderLines = ();
  $count = 0;
}

sub options {
  my ($exclude , $type , $repeatThreshold , $splitMode , $randomizations , $shuffleBed , $shuffleBed_gap ,$millionMappedReads , $tmp_split_name , $exonGtfSplits , $cluster , $species , $stranded , $bamtools , $bamdir ,   $experimentName , $pipelineDirName);
  my $spyStranded                     = 1;
  my $spySpecies                      = 1;
  my $spyBamdir                       = 1;
  my $spyBamtools                     = 1;
  my $spyExperiment                   = 1;
  my $spyPipelineDir                  = 1;
  my $spyCluster                      = 1;
  my $spyExonGtfSplits                = 1;
  my $spyTmp_split_name               = 1;
  my $spyMillionMappedReads           = 1;
  my $spyShuffleBed                   = 1;
  my $spyShuffleBed_gap               = 1;
  my $spyRandomizations               = 1;
  my $spySplitMode                    = 1;
  my $spyRepeatThreshold              = 1;
  my $spyType                         = 1;
  my $spyExclude                      = 1;



  foreach my $field (0..$#ARGV){
    if ($ARGV[$field] eq '-exclude'){
	$exclude  = $ARGV[1+$field];
	$spyExclude = 2;
	next;
    }
    if ($spyExclude == 2){
	$spyExclude = 3;
	next;
    }
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
   if ($ARGV[$field] eq '-splitMode'){
	$splitMode = $ARGV[1+$field];
	$spySplitMode = 2;
	next;
    }
    if ($spySplitMode == 2){
	$spySplitMode = 3;
	next;
    }
    if ($ARGV[$field] eq '-randomizations'){
	$randomizations = $ARGV[1+$field];
	$spyRandomizations = 2;
	next;
    }
    if ($spyRandomizations == 2){
	$spyRandomizations = 3;
	next;
    }
    if ($ARGV[$field] eq '-shuffleBed_gap'){
	$shuffleBed_gap = $ARGV[1+$field];
	$spyShuffleBed_gap = 2;
	next;
    }
    if ($spyShuffleBed_gap == 2){
	$spyShuffleBed_gap = 3;
	next;
    }
    if ($ARGV[$field] eq '-shuffleBed'){
	$shuffleBed = $ARGV[1+$field];
	$spyShuffleBed = 2;
	next;
    }
    if ($spyShuffleBed == 2){
	$spyShuffleBed = 3;
	next;
    }
    if ($ARGV[$field] eq '-millionMappedReads'){
	$millionMappedReads = $ARGV[1+$field];
	$spyMillionMappedReads = 2;
	next;
    }
    if ($spyMillionMappedReads == 2){
	$spyMillionMappedReads = 3;
	next;
    }
    if ($ARGV[$field] eq '-tmp_split_name'){
	$tmp_split_name = $ARGV[1+$field];
	$spyTmp_split_name = 2;
	next;
    }
    if ($spyTmp_split_name == 2){
	$spyTmp_split_name = 3;
	next;
    }
    if ($ARGV[$field] eq '-exonGtfSplits'){
	$exonGtfSplits = $ARGV[1+$field];
	$spyExonGtfSplits = 2;
	next;
    }
    if ($spyExonGtfSplits == 2){
	$spyExonGtfSplits = 3;
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
    if ($ARGV[$field] eq '-stranded'){
	$stranded = $ARGV[1+$field];
	$spyStranded = 2;
	next;
    }
    if ($spyStranded == 2){
	$spyStranded = 3;
	next;
    }
    if ($ARGV[$field] eq '-bamdir'){
	$bamdir = $ARGV[1+$field];
	$spyBamdir = 2;
	next;
    }
    if ($spyBamdir == 2){
	$spyBamdir = 3;
	next;
    }
    if ($ARGV[$field] eq '-bamtools'){
	$bamtools = $ARGV[1+$field];
	$spyBamtools = 2;
	next;
    }
    if ($spyBamtools == 2){
	$spyBamtools = 3;
	next;
    }
    if ($ARGV[$field] eq '-species'){
	$species = $ARGV[1+$field];
	$spySpecies = 2;
	next;
    }
    if ($spySpecies == 2){
	$spySpecies = 3;
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
  die "Error[read_statistics.pl]! Please provide the -experiment parameter\n"                             if (! defined $experimentName);
  die "Error[read_statistics.pl]! Please provide the -pipeline_dir parameter\n"                           if (! defined $pipelineDirName);
  die "Error[read_statistics.pl]! Please provide the -species parameter\n"                                if (! defined $species);
  die "Error[read_statistics.pl]! Please provide the -bamdir directory\n"                                 if (! defined $bamdir);
  die "Error[read_statistics.pl]! Please provide the -repeatThreshold parameter\n"                        if (! defined $repeatThreshold);
  die "Error[read_statistics.pl]! Please provide the -type parameter with \"bh\" (best hit) or \"all\"\n" if (! defined $type);


  #defaults
  $randomizations  = 0            unless (defined $randomizations) ;
  $stranded        = 0            unless (defined $stranded);
  $bamtools        = 'bamtools'   unless (defined $bamtools);
  $cluster         = 'off'        unless (defined $cluster);
  $exonGtfSplits   = 80           unless (defined $exonGtfSplits);
  $shuffleBed      = 'shuffleBed' unless (defined $shuffleBed);
  $splitMode       = 'off'        unless (defined $splitMode);

  return ($exclude , $type , $repeatThreshold , $splitMode , $randomizations , $shuffleBed , $shuffleBed_gap , $millionMappedReads , $tmp_split_name , $exonGtfSplits , $cluster , $species , $stranded , $bamtools , $bamdir ,   $experimentName , $pipelineDirName);
}

sub acceptedVariableSpace {
  my %space = ('-exclude' => 1 , '-type' => 1 , '-repeatThreshold' => 1 , '-splitMode' => 1 , '-randomizations' => 1 , '-shuffleBed' => 1,  '-shuffleBed_gap' => 1 , '-millionMappedReads'=> 1 , '-stranded' => 1 , '-tmp_split_name' => 1 , '-exonGtfSplits' => 1 , '-cluster' => 1 , '-experiment' => 1 , '-pipeline_dir' => 1 , '-species' => 1 , '-bamtools' => 1 , '-bamdir' => 1);
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[read_statistics.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}

sub fileNameGenerator{
  my ($nameRoot) = @_;
  my $tmp_name_counter = 0;
  my $tmp_name;
  while (!$tmp_name || -f $tmp_name) {
    $tmp_name_counter++;
    $tmp_name = "$nameRoot".".$tmp_name_counter";
  }
  return $tmp_name;
}

sub printLog {
  my ($outDir) = @_;
  open (L,">>${outDir}/log.txt") or die "Error[startPipeline.pl]! cannot create the log ${outDir}/log.txt $!\n";
  print L "#Execution started ";
  my $date = `date`;
  print L "$date";
  print L "#commandline: $0 ";
  foreach my $field (@ARGV){
    print L "$field ";
  }
  print L "\n";
  close L;
}

sub checkIndices {
  print STDERR "#check bam indices...\n";
  opendir (BD, $bamdir) or die "Error[reads_statistics.pl]! cannot  read $bamdir\n$!\n";
  while (my $file = readdir(BD)){
    next if (($file eq '.') or ($file eq '..') or ($file =~ /\.bai$/) or ($file !~ /\.bam$/));
    unless (-f "${bamdir}/${file}.bai"){
      my $cmd = "$bamtools index -in ${bamdir}/${file}";
      (system "$cmd") == 0 or die "Error[reads_statistics.pl]! cannot execute $cmd\n$!\n";
    }
  }
  closedir BD;
}

sub countAllMappingReads {
    print STDERR "#count the millions of mapping reads..\n";
    my $tot_c = 0;
    opendir (BD, $bamdir) or die "Error[reads_statistics.pl]! cannot  read $bamdir\n$!\n";
    while (my $file = readdir(BD)){
	next if (($file eq '.') or ($file eq '..') or ($file =~ /\.bai$/) or ($file !~ /\.bam$/));
	my $cmd = "$bamtools count -in $bamdir/${file}";
	my $c = `$cmd`;
	if($?){die "Error[reads_statistics.pl]! cannot execute $cmd\n$!\n";}
	chomp $c;
	$tot_c += $c;
	open (L,">>${outDir}/log.txt") or die "Error[startPipeline.pl]! cannot create the log ${outDir}/log.txt $!\n";
	print L "#file $file has $c reads\n";
	close L;
    }
    closedir (BD);
    my $millionMappedReads = $tot_c / 1000000;
    open (L,">>${outDir}/log.txt") or die "Error[startPipeline.pl]! cannot create the log ${outDir}/log.txt $!\n";
    print L "#million mapped reads $millionMappedReads\n";
    close L;
    return $millionMappedReads;
}

sub syntax_check {
  my $check1 = "$bamtools filter -region \'chr1:1..10\'";
  my $check2 = "$bamtools filter -region \'1:1..10\'";
  my $test;
  opendir (BD, $bamdir) or die "Error[reads_statistics.pl]! cannot read $bamdir\n$!\n";
  while (my $file = readdir(BD)){
    next if (($file eq '.') or ($file eq '..') or ($file =~ /\.bai$/) or ($file !~ /\.bam$/));
    $check1 .= " -in ${bamdir}/${file} > /dev/null 2> /dev/null ; echo \$?";
    $check2 .= " -in ${bamdir}/${file} > /dev/null 2> /dev/null ; echo \$?";
    $test = `$check1`;
    if ($test == 0){
      return 'chr';
    }
    $test = `$check2`;
    if ($test == 0){
      return 'noChr';
    }
    die "Error[reads_statistics.pl]! neither the command \n$check1\n nor the command\n$check2\n work.\n check manually that bamtools is working\n";
    }
}

sub set_gtf_map {
  $gtf_map = "${pipelineDirName}/experiments/${experimentName}/results/${species}.${type}.rep${repeatThreshold}".'.ex.gtf';
  if (defined $exclude){
    my $name = fileNameGenerator ('gtf_map_exclude_tmp_');
    open (F, ">$name")    or die "Error[reads_statistics.pl]! cannot generate $name $!\n";
    open (G, "<$gtf_map") or die "Error[reads_statistics.pl]! cannot read $gtf_map $!\n";
    while (my $l = <G>){
      if ($l !~/$exclude/){print F "$l"; }
    }
    close F;
    close G;
    $gtf_map = $name;
  }
  return $gtf_map;
}

sub countReadsPerExon {
  print STDERR "#count reads per exon..\n";
  open (E,"<$gtf_map") or die "Error[reads_statistics.pl]! cannot  read the $gtf_map file \n$!\n";
  while (my $line = <E>){
    if ($line=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+/){
      my $chr     = $1;
      my $start   = $2;
      my $end     = $3;
      #check chr syntax
      if (($syntaxCheck eq 'noChr') && ($chr =~/^chr/)){
	$chr=~s/^chr//;
      }
      if (($syntaxCheck eq 'chr') && ($chr !~/^chr/)){
	$chr = "chr$chr";
      }
      my $region  = "${chr}:${start}..${end}";
      my $ex_size = abs($end - $start) + 1;  # add +1 since the annotations are 1 based
      push (@allExons,$region);
      my $cmd = "$bamtools filter -region \'$region\'  ";
      $cmd .= "-isReverseStrand false " if ($stranded == 1);
      my $cmd_template = $cmd;

      my $transcript_id;
      if ($line=~/hitName \"([^\"]+)\"/){
	$transcript_id = $1;
      }
      else{
	die "Error[reads_statistics.pl]! Impossible to parse $line \n in $gtf_map \n for transcript_id\n";
      }
      #initialize
      $tx_info{$transcript_id}{'count'} = 0 unless (defined $tx_info{$transcript_id}{'count'});
      $ex_info{$region}{'count'}        = 0 unless (defined $ex_info{$region}{'count'});
      $tx_info{$transcript_id}{'size'}  = 0 unless (defined $tx_info{$transcript_id}{'size'});
      $ex_info{$region}{'size'}         = 0 unless (defined $ex_info{$region}{'size'});


      #take size
      $ex_info{$region}{'size'}         = $ex_size;
      $tx_info{$transcript_id}{'size'} += $ex_size;
      #take count
      opendir (BD, $bamdir) or die "Error[reads_statistics.pl]! cannot read $bamdir\n$!\n";
      while (my $file = readdir(BD)){
	next if (($file eq '.') or ($file eq '..') or ($file =~ /\.bai$/) or ($file !~ /\.bam$/));
	$cmd .= "-in ${bamdir}/${file} | bamtools count ";
	$count = `$cmd`; if ($?){die "Error[reads_statistics.pl]! cannot execute $cmd\n$!\n";}   chomp $count;
	$count = 0 if ($count eq "");
	$tx_info{$transcript_id}{'count'} += $count;
	$ex_info{$region}{'count'}        += $count;
	$count = 0;
	$cmd = $cmd_template;
      }
      closedir (BD);

    }
    else {die "Error[reads_statistics.pl]! Unable to parse gtf line $line \n";}
  }
  close E;
}

sub printOut {
    my ($tag) = @_;
    my $exonRpkm = fileNameGenerator("${outDir}/exons.rpkm_${tag}_");
    open (ER , ">$exonRpkm") or die "Error[reads_statistics.pl]! Cannot create $exonRpkm\n$!\n"; if ($?){exit;}
    foreach my $ex (@allExons){
	my $rpkm = $ex_info{$ex}{'count'} / ($millionMappedReads * ($ex_info{$ex}{'size'} / 1000));
	print ER "$ex" . "\t$rpkm\n";
    }
    close ER;

    my $exonRC = fileNameGenerator("${outDir}/exons.readCount_${tag}_");
    open (ER , ">$exonRC") or die "Error[reads_statistics.pl]!  Cannot create $exonRC\n$!\n"; if ($?){exit;}
    foreach my $ex (@allExons){
	my $RC = $ex_info{$ex}{'count'};
	print ER "$ex" . "\t$RC\n";
    }
    close ER;


    my $transcriptRpkm = fileNameGenerator("${outDir}/transcripts.rpkm_${tag}_");
    open (TR , ">$transcriptRpkm") or die "Error[reads_statistics.pl]! Cannot create $transcriptRpkm\n$!\n"; if ($?){exit;}
    foreach my $tx (keys %tx_info){
	my $rpkm = $tx_info{$tx}{'count'} / ($millionMappedReads * ($tx_info{$tx}{'size'} / 1000));
	print TR "$tx" . "\t$rpkm\n";
    }
    close TR;

    my $transcriptRC = fileNameGenerator("${outDir}/transcripts.readCount_${tag}_");
    open (TR , ">$transcriptRC") or die "Error[reads_statistics.pl]! Cannot create $transcriptRC\n$!\n"; if ($?){exit;}
    foreach my $tx (keys %tx_info){
	my $RC = $tx_info{$tx}{'count'};
	print TR "$tx" . "\t$RC\n";
    }
    close TR;
}

sub splitExons {
  my ($splitDir) = @_;
  (system "mkdir -p $splitDir") == 0 or die "Error[reads_statistics.pl]! cannot create the $splitDir\n$!\n";
  (system "rm -rf ${splitDir}/* 2>/dev/null") == 0 or die "Error[reads_statistics.pl]! cannot clean the $splitDir\n$!\n";
  my $allExons = `wc -l $gtf_map |  cut -f1 -d \" \"`; if($?){die "Error[reads_statistics.pl]! count the exons with\nwc -l $gtf_map | cut -f1 -d \" \"\n";}
  chomp $allExons;
  my $exonsForEachSplit = $allExons / $exonGtfSplits;
  $exonsForEachSplit    = sprintf("%.0f", $exonsForEachSplit);
  $exonsForEachSplit    = 0 if ($exonsForEachSplit == 0);

  my %transcripts;
  open (G,"<$gtf_map") or die "Error[reads_statistics.pl]! Cannot read $gtf_map\n$!\n";
  while (my $l = <G>){
      chomp $l;
      if ($l=~/hitName \"([^\"]+)\"/){
	  my $tx_id = $1;
	  push (@{$transcripts{$tx_id}} , $l);
      }
      else {die "Error[reads_statistics.pl]! Error while parsing $l\n";}
  }
  close G;

  my $c = 0;
  my $tmp_split_name = fileNameGenerator("${splitDir}/split");
  push (@allSplitNames , $tmp_split_name);
  open (TF,">$tmp_split_name") or die "Error[reads_statistics.pl]! Cannot create $tmp_split_name\n$!\n";


  foreach my $tx_id (keys %transcripts){
      foreach my $ex (@{$transcripts{$tx_id}}){
	  $c++;
	  print TF "$ex\n";
      }
      if ($c > $exonsForEachSplit){
	  close TF;
	  $tmp_split_name = fileNameGenerator("${splitDir}/split");
	  open (TF,">$tmp_split_name") or die "Error[reads_statistics.pl]! Cannot create $tmp_split_name\n$!\n";
	  push (@allSplitNames , $tmp_split_name);
	  $c = 0;
      }
  }
  close TF;

  #remove possible empty splits
  my $check_size = `wc -l $tmp_split_name |cut -f1 -d \" \" | tr -d \"\\n\"`; if($?){die "Error[reads_statistics.pl]! Cannot estimate estimate split size for $tmp_split_name\n$!\n"};
  if ($check_size == 0){
      (system "rm $tmp_split_name") == 0 or die "Error[reads_statistics.pl]! Cannot remove $tmp_split_name\n$!\n";
      pop(@allSplitNames);
  }

}


sub readClusterConfig {
  my $clusterConfigName = "${pipelineDirName}/experiments/${experimentName}/CONFIG/clusterConfig";
  open (CC,"<$clusterConfigName") or die "Error[reads_statistics.pl]! Cannot read $clusterConfigName \n$!\n";
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


sub waiting {
    my ($splitDir) = @_;
    my $time  = 15;
    my $done  = 0;
    my $total = scalar(@allSplitNames);
    while ($done < $total ){
        print STDERR "#Waiting the job to exit the cluster... " . $done . " out of " . $total . " already finished\n";
        system "sleep $time";
        $done = `ls ${splitDir}/*___read_statistics_done___ 2> /dev/null | wc -l`;
        chomp $done;
    }
    print STDERR "#All jobs exited the cluster\n";
    clearDone($splitDir);
    clusterErrorCheck();
}

sub clusterErrorCheck {
  my $clusterFileDir = "${pipelineDirName}/experiments/${experimentName}/CLUSTER_FILES/";
  my (@filesHavingError , $error , $errOut);
  my @errors = ('alloc' , 'ERROR' , 'error' , 'Error' , 'Killed' , 'memory');
  my $spy_err = 0;
  opendir (CD,"$clusterFileDir") or die "Error[reads_statistics.pl]! cannot open the directory $clusterFileDir\n$!\n";
  while( ( my $file = readdir(CD))){
    next if (($file eq '.') or ($file eq '..') or ($file !~/^split\.\d+\.sh\.e\S+$/));
    my $filename = "${clusterFileDir}/$file";
    foreach $error (@errors){
      $errOut = `grep -l $error $filename`;
      if ($errOut) {push (@filesHavingError , $errOut );$errOut = ''; $spy_err = 1; last;}
      $errOut = '';
    }
    if ($spy_err == 0){
      (system "cat $filename >> ${pipelineDirName}/experiments/${experimentName}/STDERR/bamtools") == 0 or die "Error[read_statistics.pl]! concatenate bamtools error message for $filename in STDERR $!\n";
      (system "rm $filename") == 0 or die "Error[reads_statistics.pl]! cannot remove $filename\n$!\n";
    }
    $spy_err = 0;
  }
  closedir CD;
  if (@filesHavingError){
    print "Error[reads_statistics.pl]! The following files contain the word \"error\" or \"malloc\"  or \"Killed\" or \"memory\" while running on the cluster:\n";
    foreach my $f (@filesHavingError){
      print $f ;
    }
    die "Please rerun manually the individual shell scripts before listed\n";
  }
}
sub clearDone {
    my ($splitDir) = @_;
    my $check_done_ =  `ls ${splitDir}/*___read_statistics_done___ 2> /dev/null | wc -l`;
    chomp $check_done_;
    if ($check_done_ > 0){
        system "rm ${splitDir}/*___read_statistics_done___";
    }
}

sub runOnTheCluster {
  my ($splitDir) = @_;
  print STDERR "#run on the cluster..\n";
  #create cluster scripts for each split
  readClusterConfig();
  foreach my $split (@allSplitNames){
    my $tmp_split_name = $split;
    my $spyCluster        = 1;
    my $spyRandomizations = 1;
    my $spyMillionMappedReads = 1;
    my $shellCmdLine   = "perl " . $0 . " ";
    foreach my $f (@ARGV){
    if ($f eq '-millionMappedReads'){
        $spyMillionMappedReads = 2;
        next;
    }
    if ( $spyMillionMappedReads== 2){
        $spyMillionMappedReads = 3;
        next;
    }
    if ($f eq '-cluster'){
        $spyCluster = 2;
        next;
    }
    if ($spyCluster == 2){
        $spyCluster = 3;
        next;
    }
    if ($f eq '-randomizations'){
        $spyRandomizations = 2;
        next;
    }
    if ($spyRandomizations == 2){
        $spyRandomizations = 3;
        next;
    }
    $shellCmdLine .= " $f";
    }
    $shellCmdLine .= ' -tmp_split_name  ' . "$tmp_split_name";
    $shellCmdLine .= ' -splitMode on';
    $shellCmdLine .= ' -millionMappedReads '. "$millionMappedReads";

    #create the script
    my $scriptName = "${split}.sh";
    open (S,">$scriptName") or die "Error[reads_statistics.pl]! Cannot create the cluster script for $scriptName  $!\n";
    foreach my $h (@clusterHeaderLines){
        print S "$h\n";
    }
    print S $shellCmdLine . "\n\n";
    print S "touch $scriptName"."___read_statistics_done___" . "\n";
    close S;

    #do the cmd
    my $cmd = $clusterConfigLine;
    $cmd =~s/##SCRIPT##/$scriptName/;
    system "$cmd";  if ($?){die "Error[reads_statistics.pl]! Error message returned with:\n$cmd\n$!";}
  }
}

sub concatenate {
    opendir (OD , $outDir) or die "Error[reads_statistics.pl]! cannot open $outDir\n";
    while (my $f = readdir(OD)){
	if ($f=~/^exons\.rpkm/){
	    (system "cat ${outDir}/$f >> ${outDir}/TOTexons") == 0 or die "Error[reads_statistics.pl]! cannot concatenate ${outDir}/$f\n$!\n";
	    (system "rm ${outDir}/$f") == 0 or die "Error[reads_statistics.pl]! cannot remove ${outDir}/$f\n$!\n";
	}
	if ($f=~/^transcripts\.rpkm/){
	    (system "cat ${outDir}/$f >> ${outDir}/TOTtranscripts") == 0  or die "Error[reads_statistics.pl]! cannot concatenate ${outDir}/$f\n$!\n";
	    (system "rm ${outDir}/$f") == 0 or die "Error[reads_statistics.pl]! cannot remove ${outDir}/$f\n$!\n";
	}
	if ($f=~/^exons\.readCount/){
	    (system "cat ${outDir}/$f >> ${outDir}/TOTexonsRC") == 0 or die "Error[reads_statistics.pl]! cannot concatenate ${outDir}/$f\n$!\n";
	    (system "rm ${outDir}/$f") == 0 or die "Error[reads_statistics.pl]! cannot remove ${outDir}/$f\n$!\n";
	}
	if ($f=~/^transcripts\.readCount/){
	    (system "cat ${outDir}/$f >> ${outDir}/TOTtranscriptsRC") == 0  or die "Error[reads_statistics.pl]! cannot concatenate ${outDir}/$f\n$!\n";
	    (system "rm ${outDir}/$f") == 0 or die "Error[reads_statistics.pl]! cannot remove ${outDir}/$f\n$!\n";
	}
    }
    closedir OD;
    (system "mv ${outDir}/TOTexons ${outDir}/exons.rpkm") == 0       or die "Error[reads_statistics.pl]! cannot move TOTexons\n$!\n";
    (system "mv ${outDir}/TOTtranscripts ${outDir}/transcripts.rpkm") == 0 or die "Error[reads_statistics.pl]! cannot move TOTtranscripts\n$!\n";
    (system "mv ${outDir}/TOTexonsRC ${outDir}/exons.readCount") == 0       or die "Error[reads_statistics.pl]! cannot move TOTexonsRC\n$!\n";
    (system "mv ${outDir}/TOTtranscriptsRC ${outDir}/transcripts.readCount") == 0 or die "Error[reads_statistics.pl]! cannot move TOTtranscriptsRC\n$!\n";
}

sub renaming {
    opendir (OD , $outDir) or die "Error[reads_statistics.pl]! cannot open $outDir\n$!\n";
    while (my $f = readdir(OD)){
	if ($f =~/^exons\.rpkm.+/){
	    (system "mv ${outDir}/$f $outDir/exons.rpkm") == 0 or die "Error[reads_statistics.pl]! cannot rename file ${outDir}/$f\n$!\n" ;
	}
	if ($f =~/^transcripts\.rpkm.+/){
	    (system "mv ${outDir}/$f $outDir/transcripts.rpkm") == 0 or die "Error[reads_statistics.pl]! cannot rename file ${outDir}/$f\n$!\n";
	}
	if ($f =~/^exons\.readCount.+/){
	    (system "mv ${outDir}/$f $outDir/exons.readCount") == 0 or die "Error[reads_statistics.pl]! cannot rename file ${outDir}/$f\n$!\n" ;
	}
	if ($f =~/^transcripts\.readCount.+/){
	    (system "mv ${outDir}/$f $outDir/transcripts.readCount") == 0 or die "Error[reads_statistics.pl]! cannot rename file ${outDir}/$f\n$!\n";
	}
    }
    closedir OD;
}

sub takeAverage {
    my $num = 0;
    my $tot = 0;
    open (E,"<${outDir}/exons.rpkm ") or die "Error[reads_statistics.pl]! cannot open ${outDir}/exons.rpkm\n$!\n";
    while (my $l = <E>){
	chomp $l;
	if ($l =~/^\S+\s+(\S+)/){
	    $tot += $1;
	    $num++;
	}
	else {
	    die "Error[reads_statistics.pl]! Error when trying to parse $l\n$!\n";
	}
    }
    close E;
    my $ex_avg = $tot / $num;

    $num = 0;
    $tot = 0;
    open (E,"<${outDir}/exons.readCount ") or die "Error[reads_statistics.pl]! cannot open ${outDir}/exons.readCount\n$!\n";
    while (my $l = <E>){
      chomp $l;
      if ($l =~/^\S+\s+(\S+)/){
	$tot += $1;
	$num++;
      }
      else {
	die "Error[reads_statistics.pl]! Error when trying to parse $l\n$!\n";
      }
    }
    close E;
    my $exRC_avg = $tot / $num;


    $num = 0;
    $tot = 0;
    open (T,"<${outDir}/transcripts.rpkm ") or die "Error[reads_statistics.pl]! cannot open ${outDir}/transcripts.rpkm\n$!\n";
    while (my $l = <T>){
	chomp $l;
	if ($l =~/^\S+\s+(\S+)/){
	    $tot += $1;
	    $num++;
	}
	else {
	    die "Error[reads_statistics.pl]! cannot parse $l\n$!\n";
	}
    }
    close T;
    my $tx_avg = $tot / $num;

    $num = 0;
    $tot = 0;
    open (T,"<${outDir}/transcripts.readCount ") or die "Error[reads_statistics.pl]! cannot open ${outDir}/transcripts.readCount\n$!\n";
    while (my $l = <T>){
	chomp $l;
	if ($l =~/^\S+\s+(\S+)/){
	    $tot += $1;
	    $num++;
	}
	else {
	    die "Error[reads_statistics.pl]! cannot parse $l\n$!\n";
	}
    }
    close T;
    my $txRC_avg = $tot / $num;

    open(A,">>${outDir}/average.ex.rpkm") or die "Error[reads_statistics.pl]! cannot open ${outDir}/average.ex.rpkm\n$!\n";
    print A "$ex_avg\n";
    close A;
    open(B,">>${outDir}/average.tx.rpkm") or die "Error[reads_statistics.pl]! cannot open ${outDir}/average.tx.rpkm\n$!\n";
    print B "$tx_avg\n";
    close B;
    open(A,">>${outDir}/average.ex.readCount") or die "Error[reads_statistics.pl]! cannot open ${outDir}/average.ex.readCount\n$!\n";
    print A "$exRC_avg\n";
    close A;
    open(B,">>${outDir}/average.tx.readCount") or die "Error[reads_statistics.pl]! cannot open ${outDir}/average.tx.readCount\n$!\n";
    print B "$txRC_avg\n";
    close B;
}


sub makeShufflebedChrSize {
  my (%infoAssembly) = @_;
  my $shuffleBed_chrSize = fileNameGenerator('shuffleBed_chrSize');
  open (S,">$shuffleBed_chrSize") or die "Error[reads_statistics.pl]! Cannot create $shuffleBed_chrSize $!\n";
  foreach my $chr (keys %infoAssembly){
    print S "$chr\t" .  $infoAssembly{$chr} . "\n";
  }
  close S;
  return $shuffleBed_chrSize;
}


sub takeAssemblyInfo {
  my ($assembly_dir) = @_;
  print STDERR "#take assembly info..\n";
  opendir ( DIR, $assembly_dir ) || die "Error[reads_statistics.pl]! Cannot open dir\n $assembly_dir\n";
  my %info;
  while (my $file = readdir(DIR)){
    next if (($file eq '.') or ($file eq '..'));
    open (F,"<${assembly_dir}/$file") or die "Error[reads_statistics.pl]! cannot read ${assembly_dir}/$file\n";
    my $name;
    my $totLen = 0;
    while (my $line = <F>){
      chomp $line;
      next if ($line=~/^\s*$/);
      $line =~ s/ //g;

      if ($line=~/>(\S+)/){
        $name = $1;
        next;
      }
      $totLen += length($line);
    }
    close F;
    $info{$name}  = $totLen;
  }
  closedir(DIR);
  return %info;
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
            if ($group =~/hitName \"([^\"]+)\"/){
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
  $numberTranscripts = scalar keys %infoGtfTranscript;
  return %infoGtfTranscript;
}

sub randomSortTxId {
    my (%hash) = @_;
    my (%out , @sort);
    foreach my $id (keys %hash){
	my $rNum = rand($numberTranscripts);
	$out{$rNum} = $id;
    }
    foreach my $i (sort { $b <=> $a } keys %out){
	push (@sort,$out{$i})
    }
    return @sort;
}

sub transcriptGtf2exonGtf {
    my ($in_tx_file , $infoExons_ref) = @_;
    my %ex_hash = %{$infoExons_ref};
    my %out;
    open (IN,"<$in_tx_file") or die "Error[reads_statistics.pl]! cannot read $in_tx_file\n$!\n";
    while (my $line = <IN>){
	my ($chr , $start , $strand , $tx_id);
	if ($line =~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+/){
	    $chr    = $1;
	    $start  = $2;
	    $strand = $3;
	}
	if ($line =~/hitName "([^\"]+)"/ ){
	    $tx_id = $1;
	}

	my $delta = -5;
	my ($newStart , $newEnd);
	foreach my $ex (sort {$a->{"start"} <=> $b->{"start"}}  @{$ex_hash{$tx_id}}){
	    $delta = abs($ex->{"start"} - $start) if ($delta == -5);
	    if ($start > $ex->{"start"}){
		$newStart = $ex->{"start"} + $delta;
		$newEnd   = $ex->{"end"}   + $delta;
	    }
	    if ($start <= $ex->{"start"}){
		$newStart = $ex->{"start"} - $delta;
		$newEnd   = $ex->{"end"}   - $delta;
	    }
	    my $newExon = "$chr\t" . $ex->{"source"} . "\t" . $ex->{"feature"} . "\t$newStart\t$newEnd\t" . $ex->{"score"} . "\t$strand\t" . $ex->{"frame"} . "\t" . $ex->{"group"};
	    push (@{$out{$tx_id}},$newExon);
	}
    }
    close IN;
    return %out;
}








sub takeBlocks {
    my ($current_gtf) = @_;
    my $cmd_sort = 'sort -k 1,1 -k 7,7 -k 4,4n  ' . "$current_gtf";
    my @sorted_gtf = `$cmd_sort`; if ($?) {die "Error[reads_statistics.pl]! cannot run $cmd_sort\n";}
    my ($oldChr , $oldStart , $oldEnd , $oldStrand , %extractedTX , $foundExons , $exonsToBeFound , @current_block  );
    #initialize
    my $old_l = shift(@sorted_gtf);
    my $c = 0;
    if ($old_l =~/hitName \"([^\"]+)\"/){ 
	my $tx_id = $1;
	$extractedTX{$tx_id} = 1;
	push (@current_block , $tx_id);
	$exonsToBeFound += scalar(@{$gtf_map_info_ex{$tx_id}});
	$foundExons++;
    }
    my $overlap = 0;

    #LOOP
    while (my $l = shift(@sorted_gtf)){
	my $tx_id;
	if ($l =~/hitName \"([^\"]+)\"/){
	    $tx_id = $1;
	}

	if ((checkExonOverlap($l,$old_l) == 1) or ($foundExons < $exonsToBeFound)){
	    if (! defined $extractedTX{$tx_id}){ 
		$extractedTX{$tx_id} = 1;
		push (@current_block , $tx_id);
		$exonsToBeFound += scalar(@{$gtf_map_info_ex{$tx_id}});
	    }
	    $foundExons++;
	    $old_l      = $l;
	    next;
	}
	#save and reinitialize
	push (@{$blocks[$c]} , @current_block );
	$c++;
	@current_block             = ();
	%extractedTX               = ();
	$exonsToBeFound            = 0;
	push (@current_block , $tx_id);
	$exonsToBeFound += scalar(@{$gtf_map_info_ex{$tx_id}});
	$extractedTX{$tx_id} = 1;
	$foundExons          = 1;
	$old_l               = $l;
    }
    push (@{$blocks[$c]} , @current_block );
}
sub blockGtf2exonGtf {
    my ($blockFile , $outFile) = @_;
    open (B,"<$blockFile") or die "Error[reads_statistics.pl]! cannot open $blockFile $!\n";
    open (O,">$outFile") or die "Error[reads_statistics.pl]! cannot open $outFile $!\n";

    while (my $current_block = <B>){
	#take shuffled block coordinates
	my ($p_chr , $p_start , $p_strand , $block_in);
	if ($current_block =~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+block \"([^\"]+)\";/){
	    $p_chr    = $1;
	    $p_start  = $2;
	    $p_strand = $3;
	    $block_in = $4;
	}

	#to do
	my $projection_delta = $p_start - $blockStarts[$block_in];

	#take internal transcript deltas
	foreach my $tx_id (@{$blocks[$block_in]}){
	    foreach my $ex_index (0..$#{$deltas4block[$block_in]{$tx_id}}){
		my $new_chr   = $p_chr;
		my $source    = $gtf_map_info_ex{$tx_id}[$ex_index]->{'source'};
		my $feature   = 'exon';
		my $new_start = $deltas4block[$block_in]{$tx_id}[$ex_index]{'delta_start'} + $p_start;
		my $new_end   = $deltas4block[$block_in]{$tx_id}[$ex_index]{'delta_end'}   + $p_start;
		my $score     = $gtf_map_info_ex{$tx_id}[$ex_index]->{'score'};
		my $strand    = $p_strand;
		my $frame     = $gtf_map_info_ex{$tx_id}[$ex_index]->{'frame'};
		my $group     = $gtf_map_info_ex{$tx_id}[$ex_index]->{'group'};
		print O "$new_chr\t$source\t$feature\t$new_start\t$new_end\t$score\t$strand\t$frame\t$group\n";
	    }

	}

    }
    close B;
    close O;
}

sub doBlocks {
    my $name = fileNameGenerator ('map.block.gtf');
    open (O,">$name") or die "Error[reads_statistics]! cannot create $name\n$!\n";
    foreach my $b (0..$#blocks){
	#print "block $b ";
	my (@allCoordinates , $strand , $chr , $outLine);
	foreach my $tx_id (@{$blocks[$b]}){
	    push (@allCoordinates , $gtf_map_info_tx{$tx_id}->{'start'});
	    push (@allCoordinates , $gtf_map_info_tx{$tx_id}->{'end'});
	    $strand = $gtf_map_info_tx{$tx_id}->{'strand'};
	    $chr    = $gtf_map_info_tx{$tx_id}->{'chr'};
	}
	my $blockStart  = min (@allCoordinates);
	my $blockEnd    = max (@allCoordinates);

	$outLine .= $chr              . "\t";
	$outLine .= 'read_statistics' . "\t";
	$outLine .= "block"           . "\t";
	$outLine .= $blockStart       . "\t";
	$outLine .= $blockEnd         . "\t";
	$outLine .= '.'               . "\t";
	$outLine .= $strand           . "\t";
	$outLine .= '.'               . "\t";
	$outLine .= "block \"$b\";";
	print O "$outLine\n";

	#take block internal deltas
	foreach my $tx_id (@{$blocks[$b]}){
	    foreach my $ex (@{$gtf_map_info_ex{$tx_id}}){
		my %hash;
		$hash{'delta_start'} = $ex->{'start'} -  $blockStart;
		$hash{'delta_end'}   = $ex->{'end'} -  $blockStart;
		push (@{$deltas4block[$b]{$tx_id}} , \%hash);
	    }
	}

	#save block start
	$blockStarts[$b] = $blockStart;
    }
    close O;
    return $name;

}
sub checkExonOverlap {
    my ($old , $new) = @_;
    my ($o_chr , $o_start , $o_end , $o_strand , $n_chr , $n_start , $n_end , $n_strand);
    if ($old=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+/){
	$o_chr    = $1;
	$o_start  = $2;
	$o_end    = $3;
	$o_strand = $4;
    }
    if ($new=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+/){
	$n_chr    = $1;
	$n_start  = $2;
	$n_end    = $3;
	$n_strand = $4;
    }
    my $overlap = 0;
    if (($o_chr eq $n_chr) && ($o_strand eq $n_strand)){
	$overlap = verifyOverlap ($o_start , $o_end , $n_start , $n_end);
    }

    return $overlap;
}
sub verifyOverlap {
  my ($traSta ,  $traEnd , $geneSta , $geneEnd ) = @_;
  my $overlap = 0;
  if ((($traSta <= $geneEnd) && ($traSta >= $geneSta)) || (($traEnd >= $geneSta) && ($traEnd <= $geneEnd)) || (($traSta <= $geneSta) && ($traEnd >= $geneEnd)) ||   (($traSta >= $geneSta) && ($traEnd <= $geneEnd))){
    $overlap = 1;
  }
  return $overlap;
}
sub createRandomWithShuffleBed {
  my ($ran)       = @_;
  my $tmp = fileNameGenerator("_tmp_");
  my $cmd = "$shuffleBed -i $gtf_map_block -g $shuffleBed_chrSize";
  $cmd .= " -excl $shuffleBed_gap " if (defined $shuffleBed_gap);
  $cmd .= " > $tmp" ;
  (system "$cmd") == 0 or die "Error[reads_statistics.pl]! cannot run $cmd \n$!\n";            # system "cat $tmp";
  blockGtf2exonGtf($tmp , $ran);
  system "rm $tmp";
}















######old

#	#take block
#	while (($foundExons < $exonsToBeFound) or ( $overlap == 1 )){
#	    $foundExons++ if ($foundExons == -1);
#	    $overlap = checkExonOverlap($l,$old_l);
#	    if ($l =~/transcript_id \"([^\"]+)\"/){
#		my $tx_id = $1;
#		if (! defined $extractedTX{$tx_id}){ 
#		    $extractedTX{$tx_id} = 1;
#		    push (@current_block , $tx_id);
#		    $exonsToBeFound += scalar(@{$gtf_map_info_ex{$tx_id}});
#		}	      
#		$foundExons++;
#	    }
#	    $old_l = $l;
#	}

	#take delta
#	my %delta4currentBlockMembers;
#	if (scalar @current_block > 1){
#	    my $base_tx = $current_block[0];
#	    my $base_start = $gtf_map_info_tx{$base_tx}{'start'};
#	    foreach my $i (1..$#current_block){
#		my $other_tx    = $current_block[$i];
#		my $other_start = $gtf_map_info_tx{$other_tx}{'start'};
#		$delta4currentBlockMembers{$other_tx} = $other_start ;
#	    }
#	}

	#save and reinitialize
#	push (@blocks , \@current_block);
#	push (@deltas4block , \%delta4currentBlockMembers);
#	@current_block             = ();
#	%extractedTX               = ();
#	%delta4currentBlockMembers = ();
#	$exonsToBeFound = 0;
#	$foundExons     = -1;
#   }
#}


#sub createRandomWithShuffleBed2 {
#  my ($ran)       = @_;
#  my $tmp_cov     = 0;
#  my $switch      = 0;
#  while ($switch == 0){
#      my $tmp = fileNameGenerator("_tmp_");
#      my @sortedTxId  = randomSortTxId(%gtf_map_info_tx);
#      my $cmd = "$shuffleBed -i $gtf_map_block -g $shuffleBed_chrSize"; 
#      $cmd .= " -excl $shuffleBed_gap " if (defined $shuffleBed_gap); 
#      $cmd .= " > $tmp" ;
#      (system "$cmd") == 0 or die "Error[reads_statistics.pl]! cannot run $cmd \n$!\n";            # system "cat $tmp";
#      my %ran_ex = transcriptGtf2exonGtf($tmp , \%gtf_map_info_ex);
#      (system "rm $tmp") == 0 or die "Error[reads_statistics.pl]! cannot remove $tmp \n$!\n";
#      open (R,">>$ran")  or die "Error[reads_statistics.pl]! cannot removeopen $ran \n$!\n";
#      while (($tmp_cov < $real_genome_coverage)&&(my $randomTx = shift (@sortedTxId))){
#	  foreach my $e (@{$ran_ex{$randomTx}}){
#	      print R $e . "\n";
#	  }
#	                                                                                            #system "cat $ran";
#	  $tmp_cov = takeCoverage($ran);                                                            #print "$tmp_cov  $real_genome_coverage\n";
 #     }
#      close R;      
#      $switch = 1 unless ($tmp_cov < $real_genome_coverage);
#  }
#}


#sub doTranscriptGtf {
#    my (%infoTranscripts) = @_;
#    my $name = fileNameGenerator ('map.tx.gtf');
#    open (O,">$name") or die "Error[reads_statistics]! cannot create $name\n$!\n";
#    foreach my $id (keys %infoTranscripts){
#	my $outLine;
#	$outLine .= $infoTranscripts{$id}->{'chr'}     . "\t";
#	$outLine .= $infoTranscripts{$id}->{'source'}  . "\t";
#	$outLine .= "transcript" . "\t";
#	$outLine .= $infoTranscripts{$id}->{'start'}   . "\t";
#	$outLine .= $infoTranscripts{$id}->{'end'}     . "\t";
#	$outLine .= $infoTranscripts{$id}->{'score'}   . "\t";
#	$outLine .= $infoTranscripts{$id}->{'strand'}  . "\t";
#	$outLine .= $infoTranscripts{$id}->{'frame'}   . "\t";
#	$outLine .= $infoTranscripts{$id}->{'group'}   . "\t";
#	print O "$outLine\n";
#    }
#    close O;
#    return $name;
#}


#sub takeDelta {
#    my (@current_block) = @_;
#    my %hash;
#    my $base_tx = $current_block[0];
#    my $base_start = $gtf_map_info_tx{$base_tx}{'start'};
#    foreach my $i (0..$#current_block){
#	my $other_tx     = $current_block[$i];
#	my $other_start  = $gtf_map_info_tx{$other_tx}{'start'};
#	$hash{$other_tx} = $other_start - $base_start;
#    }
#    return %hash;
#}
