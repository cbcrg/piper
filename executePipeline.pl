#!/usr/bin/env perl
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
my ($stepName , $experimentName , $pipelineDirName) = options();
my $pipelineFileName = "${pipelineDirName}/experiments/${experimentName}/RNAmapping_pipeline.txt";
my @allSteps = ("pre_processing","database","blast","exonerate","prepare_mfa","alignment","similarity","csf");
my $endLine = takeEndLine();


open (P,"<$pipelineFileName") or die "Error[executePipeline.pl]! Cannot read the $pipelineFileName file $!\n";
my @pipelineFile = <P>;
close P;
my $tmpFileName = fileNameGenerator('runPipeline');
open (R,">$tmpFileName") or die "Error[executePipeline.pl]! Cannot create the file $tmpFileName\n$!\n";
my $parameters = 0;
my $workingDir;
foreach my $line (@pipelineFile){
    chomp $line;
    $parameters = 1 if ($line=~/^############$/);
    print R "$line\n" if ($parameters == 0);
    if ($line=~/^PIPELINEDIR=(\S+)$/){$workingDir = $1;}
}


my $startLine = takeStartLine();
sanityCheck($startLine , $endLine);
my $startSpy = 0;
foreach my $line (@pipelineFile){
    chomp $line;
    $startSpy = 1 if ($line=~/$startLine/);
    $startSpy = 0 if ($line=~/$endLine/);
    next if ($startSpy == 0);
    print R "$line\n";
}

close R;





#RUN IT!
#die "check!\n";
(system "sh $tmpFileName") == 0 or die "Error[executePipeline]! Cannot run $tmpFileName $!\n";
system "rm $tmpFileName";





###
#FUNCTIONS
sub takeEndLine {
    my $endLine;
    if ($stepName eq 'pre_processing'){
	$endLine = '#CREATE DATABASE';
    }
    if ($stepName eq 'database'){
	$endLine = '#BLAST SEARCH';
    }
    if ($stepName eq 'blast'){
	$endLine = '#EXONERATE REMAPPING';
    }
    if ($stepName eq 'exonerate'){
	$endLine = '#PREPARE_MFA';
    }
    if ($stepName eq 'prepare_mfa'){
	$endLine = '#MULTIPLE SEQUENCE ALIGNMENT';
    }
    if ($stepName eq 'alignment'){
	$endLine = '#EXTRACT SIMILARITY';
    }
    if ($stepName eq 'similarity'){
	$endLine = '#CSF';
    }
    if ($stepName eq 'csf'){
	$endLine = '#END';
    }
    return $endLine;
}
 sub takeStartLine {
   my $startLine = '############';
   my $stop = 0;
   foreach my $s (@allSteps){
     if ($stop == 1){last;}
     if ($s eq $stepName){$stop=1;}

      if ($s eq 'pre_processing'){
        next;
      }

      if ($s eq 'database'){
        next;
      }

     if ($s eq 'blast'){
       my $check = `ls ${workingDir}/experiments/$experimentName/BLAST_OUT | wc -l`;
       chomp $check;
       if ($check > 0) 	{
	 $startLine = '#EXONERATE REMAPPING';
	 next;
       }
       else {
	 $startLine = '############';
	 last;
       }
     }

     if ($s eq 'exonerate'){
       my $check = `ls ${workingDir}/experiments/$experimentName/EXONERATE_OUT | wc -l`;
       chomp $check;
       if ($check > 0) 	{
	 $startLine = '#PREPARE_MFA';
	 next;
       }
       else {
	 $startLine = '#EXONERATE REMAPPING';
	 last;
       }
     }



     if ($s eq 'prepare_mfa'){
       my $check = `ls ${workingDir}/experiments/$experimentName/multifasta4EachTx | wc -l`;
       chomp $check;
       if ($check > 0) 	{
	 $startLine = '#MULTIPLE SEQUENCE ALIGNMENT';
	 next;
       }
       else {
	 $startLine = '#PREPARE_MFA';
	 last;
       }
     }


     if ($s eq 'alignment'){
       my $check = `ls ${workingDir}/experiments/$experimentName/outAlignments | wc -l`;
       chomp $check;
       if ($check > 0) 	{
	 $startLine = '#EXTRACT SIMILARITY';
	 next;
       }
       else {
	 $startLine = '#MULTIPLE SEQUENCE ALIGNMENT';
	 last;
       }
     }


     if ($s eq 'similarity'){
       my $check = `ls ${workingDir}/experiments/$experimentName/outSim | wc -l`;
       chomp $check;
       if ($check > 0) 	{
	 $startLine = '#CSF';
	 next;
       }
       else {
	 $startLine = '#EXTRACT SIMILARITY';
	 last;
       }
     }

     if ($s eq 'csf'){
       my $check = `ls ${workingDir}/experiments/$experimentName/fasta_aln | wc -l`;
       chomp $check;
       if ($check > 0) 	{
	 $startLine = '#END';
	 next;
       }
       else {
	 $startLine = '#CSF';
	 last;
       }
     }

   }
   return $startLine;
 }

sub acceptedVariableSpace {
  my %space = ('-step' => 1 , '-experiment' => 1 , '-pipeline_dir' => 1);
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[executePipeline.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}

sub options {
  my ($stepName , $experimentName , $pipelineDirName);
  my $spyStep          = 1;
  my $spyExperiment    = 1;
  my $spyPipelineDir   = 1;
  foreach my $field (0..$#ARGV){
    if ($ARGV[$field] eq '-step'){
	$stepName = $ARGV[1+$field];
	$spyStep = 2;
	next;
    }
    if ($spyStep == 2){
	$spyStep = 3;
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

  die "Error[executePipeline.pl]! You must provide the -pipeline_dir parameter.\n" if ($spyPipelineDir != 3);
  die "Error[executePipeline.pl]! Accepted steps are pre_processing|database|blast|exonerate|prepare_mfa|alignment|similarity \n" if ((defined $stepName) && ($stepName ne 'pre_processing') && ($stepName ne 'database') && ($stepName ne 'blast') && ($stepName ne 'exonerate') && ($stepName ne 'prepare_mfa') && ($stepName ne 'alignment') && ($stepName ne 'similarity') && ($stepName ne 'csf'));
  $stepName = 'exonerate' if (! defined $stepName);
  if ($spyExperiment  != 3){print "Error[executePipeline.pl]! The -experiment parameter is mandatory. Please choose your experiment\n"; system "ls ${pipelineDirName}/experiments/ "; exit;}


  return ($stepName , $experimentName , $pipelineDirName);
}


sub fileNameGenerator{
  my ($nameRoot) = @_;
  my $tmp_name_counter = 0;
  my $tmp_name;
  while (!$tmp_name || -f $tmp_name) {
    $tmp_name_counter++;
    $tmp_name = "$nameRoot$$".".$tmp_name_counter";
  }
  return $tmp_name;
}

sub sanityCheck{
  my ($startLine,$endLine) = @_;
  my %steps = (
	       '############'                 => 1,
	       '#CREATE DATABASE'             => 2,
	       '#BLAST SEARCH'                => 3,
	       '#EXONERATE REMAPPING'         => 4,
	       '#PREPARE_MFA'                 => 5,
	       '#MULTIPLE SEQUENCE ALIGNMENT' => 6,
	       '#EXTRACT SIMILARITY'          => 7,
	       '#CSF'                         => 8,
	       '#END'                         => 9
	      );
  if ($steps{$startLine} >= $steps{$endLine}){close R; system "rm $tmpFileName"; die "Error[executePipeline]! It looks like that the pipeline was already run to $stepName step\n";};
  my %check = (
	       '1' => "############",
	       '2' => "${pipelineDirName}/allGenomeInfo/",
	       '3' => "${pipelineDirName}/experiments/${experimentName}/BLAST_OUT/",
	       '4' => "${pipelineDirName}/experiments/${experimentName}/EXONERATE_OUT/",
	       '5' => "${pipelineDirName}/experiments/${experimentName}/multifasta4EachTx/",
	       '6' => "${pipelineDirName}/experiments/${experimentName}/outAlignments/",
	       '7' => "${pipelineDirName}/experiments/${experimentName}/outSim/",
	       '8' => "${pipelineDirName}/experiments/${experimentName}/fasta_aln/",
	       '9' => ""
	      );
  foreach my $d ($steps{$startLine}..8){
    next if (($d == 1) or ($d == 2));
    my $dir      = $check{$d};
    my $checkRes = `ls $dir 2> /dev/null | wc -l`;
    chomp $checkRes;
    if ($checkRes > 0){
      close R;
      system "rm $tmpFileName";
      die "Error[executePipeline]! It looks like that the pipeline was already run to $check{$d} step\n";
    }
  }

}



sub help_message {
my $helpMessage = "\nNAME
executePipeline.pl - Run the pipeline on a specified experiment up to the step indicated. It skips the steps that are already been performed\n
SYNOPSIS
./executePipeline.pl -pipeline_dir -experiment [-step pre_processing|database|blast|exonerate|prepare_mfa|alignment|similarity|csf]\n
DESCRIPTION
   * executePipeline.pl run the pipeline created by startPipeline.pl
   * It runs the pipeline to the experiment and to the step selected by the user
   * This script check that the steps already done are skipped
OPTIONS
   * -step [pre_processing|database|blast|exonerate|prepare_mfa|alignment|similarity|csf] <Default: exonerate>
";
print "$helpMessage\n\n";
exit;
}






