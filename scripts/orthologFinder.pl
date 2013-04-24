#!/usr/bin/perl -w
use strict;
use warnings;
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
my ($query , $blastName , $strategyName , $xdformatName , $referenceGenomeName , $query_gtf_name , $experimentName , $pipelineDirName) = options();
my $configName         = "${pipelineDirName}/experiments/${experimentName}/CONFIG/blastConfig";

#TARGET DB
my $db_folder;
if (($strategyName eq 'wublastn') or ($strategyName eq 'wublastn_opt')){$db_folder = 'wublastn_db';}
if (($strategyName eq 'abblastn') or ($strategyName eq 'abblastn_opt')){$db_folder = 'abblastn_db';}
if ($strategyName eq 'wublastr'){$db_folder = 'wublastr_db';}
if ($strategyName eq 'abblastr'){$db_folder = 'abblastr_db';}

#CHECK
die "Error[orthologFinder.pl]! cannot open $pipelineDirName/experiments/$experimentName/pre_processing/inputInfo/allGenomeInfo/\nTo run orthologFinder you need this folder. In case the reference species is in allGenomeInfo, you cannot consider it anyway because it could be formatted in a strange way (i.e. splitByChromosome). So you always need the preprocessing folder $!\n" unless (-d "$pipelineDirName/experiments/$experimentName/pre_processing/inputInfo/allGenomeInfo/");
#BLAST
my $referenceName = `ls $pipelineDirName/experiments/$experimentName/pre_processing/inputInfo/allGenomeInfo/`;
chomp $referenceName;
my $targetDB = "$pipelineDirName/experiments/$experimentName/pre_processing/inputInfo/allGenomeInfo/$referenceName/${db_folder}/db";

my $configLine;
readConfig();
my $blast_cmd = "$blastName $targetDB $query ";
$blast_cmd .= " -mformat=2 -V=1 -B=1 -spoutmax=1 -warnings -errors -notes ";
$blast_cmd .= " $configLine"  if (defined $configLine);
#my $blastOut = `$blast_cmd | sort -gk3 | head -1`; # if ($?){ die "Error[orthologFinder.pl]! error with command:\n$blast_cmd\n$!\n"; }
my $blastOut = `$blast_cmd`;
if ($blastOut){
  chomp $blastOut;
}
else {
 print "0";
 exit;
}


#TAKE REAL QUERY COORDINATES
open (Q,"<$query") or die "Error[orthologFinder.pl]! cannot open $query\n$!\n";
my $queryName = <Q>;
close Q;
chomp $queryName;
$queryName =~ s/^>//;
my @queryName_gtfLine_allExons = `grep $queryName $query_gtf_name`; if ($?){ die "Error[orthologFinder.pl]! error with command:\ngrep $queryName $query_gtf_name\n$!\n"; }
if (! @queryName_gtfLine_allExons){
  die "Error[orthologFinder.pl]! the $query_gtf_name gtf file does not contain the query $queryName\n";
}







#COMPARE
my ($b_chr , $b_s_start , $b_s_end , $b_q_frame , $b_s_frame , $b_strand);
if ($blastOut=~/^\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s*$/){
  $b_chr     = $1;
  $b_q_frame = $2;
  $b_s_frame = $3;
  $b_s_start = $4;
  $b_s_end   = $5;

  $b_strand = '+' if (($b_q_frame=~/\+/) && ($b_s_frame=~/\+/));
  $b_strand = '-' if (($b_q_frame=~/\+/) && ($b_s_frame=~/-/));
  $b_strand = '-' if (($b_q_frame=~/-/) && ($b_s_frame=~/\+/));
  $b_strand = '+' if (($b_q_frame=~/-/) && ($b_s_frame=~/-/));
}
else {
  die "Error[orthologFinder.pl]! impossible to parse $blastOut\n";
}
foreach my $queryName_gtfLine (@queryName_gtfLine_allExons){
  my ($g_chr , $g_start , $g_end , $g_strand);
  if ($queryName_gtfLine=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)/){
    $g_chr    = $1;
    $g_start  = $2;
    $g_end    = $3;
    $g_strand = $4;
  }
  else {
    die "Error[orthologFinder.pl]! $queryName_gtfLine is a wrong gtf format\n";
  }
  if (($g_chr eq $b_chr) && ($g_strand eq $b_strand) ){
    verifyOverlap($b_s_start , $b_s_end , $g_start , $g_end );
  }
  else {
    print "0";
    exit;
  }
}

print "0";
exit;


















#FUNCTIONS
sub help_message {
my $helpMessage = "\nNAME
orthologFinder.pl - Perform a reciprocal blast to establish the ortholog relationship among the exonerated output and the original query

SYNOPSIS
orthologFinder.pl -query -experiment -pipeline_dir [-blast -blast_strategy -xdformat -reference_genome -query_gtf]

DESCRIPTION
   * orthologFinder.pl takes the exonerated output returned by exonerateRemapping.pl blast it VS the reference genome
   * the blast and the parameters will be the same to the ones used to screen the query.

   * If the reference genome wasn't yet formatted (by the optional pre_processing step) this must be provided via: -reference_genome multi-FASTAfile
   * In this case the input reference genome dill be formatted and created at: experimentName/pre_processing/inputInfo/allGenomeInfo

   * By deafault it will be used the transcript .gtf file generated during the pre_processing step at: experimentName/pre_processing/inputInfo/
   * If the pre_processing was skipped the gtf file containing the query positions on the reference genome must be provided via: -query_gtf transcript.gtf

   OPTIONS
   * -blast                                                                          <Default: wu-blastn>
   * -xdformat                                                                       <Default: xdformat>
   * -blast_strategy [wublastn|abblastn|wublastr|abblastr|wublastn_opt|abblastn_opt] <Default: wublastn_opt>
   * -reference_genome                                                               <Default: none>
   * -query_gtf                                                                      <Default: none>

TROUBLESHOOTING
  * The user cannot add to the blastConfig configuration file the parameters: -mformat -V -B -spoutmax -warnings -notes -errors
  * This is because these parameters are used by orthologFinder.pl to call blast and therefore are not editable by the user
";
print "$helpMessage\n\n\n";
exit;
}


sub options {
  my ($query , $blastName , $strategyName , $xdformatName , $referenceGenomeName , $query_gtf_name , $experimentName , $pipelineDirName);
  my $spyQuery                = 1;
  my $spyBlast                = 1;
  my $spyStrategy             = 1;
  my $spyXdformat             = 1;
  my $spyReferenceGenome      = 1;
  my $spyQuery_gtf_name       = 1;
  my $spyExperiment           = 1;
  my $spyPipelineDir          = 1;

  foreach my $field (0..$#ARGV){
   if ($ARGV[$field] eq '-query'){
      $query = $ARGV[1+$field];
      $spyQuery = 2;
      next;
    }
    if ($spyQuery == 2){
      $spyQuery = 3;
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
    if ($ARGV[$field] eq '-blast_strategy'){
	$strategyName = $ARGV[1+$field];
	$spyStrategy = 2;
	next;
    }
    if ($spyStrategy == 2){
	$spyStrategy = 3;
	next;
    }
    if ($ARGV[$field] eq '-xdformat'){
      $xdformatName = $ARGV[1+$field];
      $spyXdformat = 2;
      next;
    }
    if ($spyXdformat == 2){
      $spyXdformat = 3;
      next;
    }
   if ($ARGV[$field] eq '-reference_genome'){
	$referenceGenomeName = $ARGV[1+$field];
	$spyReferenceGenome = 2;
	next;
    }
    if ($spyReferenceGenome == 2){
	$spyReferenceGenome = 3;
	next;
    }
   if ($ARGV[$field] eq '-query_gtf'){
	$query_gtf_name = $ARGV[1+$field];
	$spyQuery_gtf_name = 2;
	next;
    }
    if ($spyQuery_gtf_name == 2){
	$spyQuery_gtf_name = 3;
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

  if ((defined $strategyName) && ($strategyName ne "wublastn") && ($strategyName ne "abblastn") && ($strategyName ne "wublastr") && ($strategyName ne "abblastr") && ($strategyName ne "wublastn_opt") && ($strategyName ne "abblastn_opt") ){die "Error[orthologFinder.pl]! The supported blast flavors are wublastn abblastn wublastr abblastr wublastn_opt abblastn_opt\n";}
  $xdformatName           = 'xdformat'                                                                          if (! defined $xdformatName);
  $query_gtf_name         = 'none'                                                                              if (! defined $query_gtf_name);
  $referenceGenomeName    = 'none'                                                                              if (! defined $referenceGenomeName);
  $blastName              = 'wu-blastn'                                                                         if (! defined $blastName);
  $strategyName           = 'wublastn_opt'                                                                      if (! defined $strategyName);
  die "Error[orthologFinder.pl]! You must provide the -query parameter. Please indicate the FASTA query file\n" if ($spyQuery != 3);
  die "Error[orthologFinder.pl]! Please provide the -experiment parameter\n"                                    if (! defined $experimentName);
  die "Error[orthologFinder.pl]! Please provide the -pipeline_dir parameter\n"                                  if (! defined $pipelineDirName);

  die "Error[orthologFinder.pl]! Neither the  ___tx___.gtf was found in the pre_processing folder, nor the query.gtf file was provided. Please specify the query gtf file via -query_gtf\n" if (($query_gtf_name eq 'none') && (! -f "$pipelineDirName/experiments/$experimentName/pre_processing/inputInfo/___tx___.gtf"));
  die "Error[orthologFinder.pl]! Neither the  reference allGenomeInfo was found in the pre_processing folder, nor the referenceGenome file was provided. Please specify the query gtf file via -reference_genome\n" if (($referenceGenomeName eq 'none') && (! -d "$pipelineDirName/experiments/$experimentName/pre_processing/inputInfo/allGenomeInfo/"));

  if (-f "$pipelineDirName/experiments/$experimentName/pre_processing/inputInfo/___tx___.gtf") {$query_gtf_name = "$pipelineDirName/experiments/$experimentName/pre_processing/inputInfo/___tx___.gtf";}
  die "Error[orthologFinder.pl]! You must provide the -reference_genome to make the reciprocal blast if you did not ran the pre_processing step before\n" if ((! -d "$pipelineDirName/experiments/$experimentName/pre_processing/inputInfo/allGenomeInfo/") && ($referenceGenomeName eq 'none'));

  return ($query , $blastName , $strategyName , $xdformatName , $referenceGenomeName , $query_gtf_name , $experimentName , $pipelineDirName);
}

sub acceptedVariableSpace {
  my %space = ('-query' => 1  , '-blast' => 1  , '-blast_strategy' => 1  , '-xdformat' => 1  , '-reference_genome' => 1  , '-query_gtf' => 1  ,  '-experiment' => 1 , '-pipeline_dir' => 1  );
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[prepare_mfa4alignment.pl]! $ARGV[$field] it is not a valid parameter\n";
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
    $tmp_name = "${nameRoot}_$$".".$tmp_name_counter";
  }
  return $tmp_name;
}

sub readConfig {
    open (C,"<$configName") or die "Error[orthologFinder.pl]! Cannot read the blast configuration file $!\n";
    foreach my $line (<C>){
	chomp $line;
	next if ($line=~/^\s*$/);
	$configLine .= $line . " ";
    }
    close C;
    die "Error[orthologFinder.pl]! $configName contains the -V parameter. Please remove it to run orthologFinder.pl\n" if (($configLine=~/-V /) or ($configLine=~/ V /) or ($configLine=~/-V=/) or ($configLine=~/V=/));
    die "Error[orthologFinder.pl]! $configName contains the -B parameter. Please remove it to run orthologFinder.pl\n" if (($configLine=~/-B /) or ($configLine=~/ B /) or ($configLine=~/-B=/) or ($configLine=~/B=/));
    die "Error[orthologFinder.pl]! $configName contains the -spoutmax parameter. Please remove it to run orthologFinder.pl\n" if (($configLine=~/-spoutmax /) or ($configLine=~/ spoutmax /) or ($configLine=~/-spoutmax=/) or ($configLine=~/spoutmax=/));
    die "Error[orthologFinder.pl]! $configName contains the -mformat parameter. Please remove it to run orthologFinder.pl\n" if (($configLine=~/-mformat /) or ($configLine=~/ mformat /) or ($configLine=~/-mformat=/) or ($configLine=~/mformat=/));
    die "Error[orthologFinder.pl]! $configName contains the -warnings parameter. Please remove it to run orthologFinder.pl\n" if (($configLine=~/-warnings/) or ($configLine=~/ warnings /));
    die "Error[orthologFinder.pl]! $configName contains the -notes parameter. Please remove it to run orthologFinder.pl\n" if (($configLine=~/-notes/) or ($configLine=~/ notes /));
    die "Error[orthologFinder.pl]! $configName contains the -errors parameter. Please remove it to run orthologFinder.pl\n" if (($configLine=~/-errors/) or ($configLine=~/ errors /));
}

sub verifyOverlap {
  my ($b_start , $b_end , $g_start , $g_end ) = @_;
  if ((($b_start <= $g_end) && ($b_start >= $g_start)) || (($b_end >= $g_start) && ($b_end <= $g_end)) || (($b_start <= $g_start) && ($b_end >= $g_end)) ||   (($b_start >= $g_start) && ($b_end <= $g_end))){
    print "1";
    exit;
  }
}
