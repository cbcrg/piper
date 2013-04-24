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
my ($tcoffeeName  , $phylocsfName , $experimentName , $pipelineDirName , $phyloCSFparametersName , $referencePhyloCsfName ) = options();
my $databaseName       = "${pipelineDirName}/allGenomeInfo/";
my $fasta_aln          = "${pipelineDirName}/experiments/${experimentName}/fasta_aln/";
my $outAlignments      = "${pipelineDirName}/experiments/${experimentName}/outAlignments/";

my %phylogenyNames;
if ($phyloCSFparametersName eq '29mammals'){
    %phylogenyNames = (
	"human"        => "Human"         ,
	"chimp"        => "Chimp"         ,

	"macaque"      => "Rhesus"        ,
	"rhesus"       => "Rhesus"        ,
	"tarsier"      => "Tarsier"       ,
	"mouse_lemur"  => "Mouse_lemur"   ,
	"bushbaby"     => "Bushbaby"      ,
	"treeshrew"    => "TreeShrew"     ,
	"mouse"        => "Mouse"         ,
	"rat"          => "Rat"           ,
	"kangaroo_rat" => "Kangaroo_rat"  ,
	"guinea_pig"   => "Guinea_Pig"    ,
	"squirrel"     => "Squirrel"      ,
	"rabbit"       => "Rabbit"        ,
	"pika"         => "Pika"          ,
	"alpaca"       => "Alpaca"        ,
	"dolphin"      => "Dolphin"       ,
	"cow"          => "Cow"           ,
	"horse"        => "Horse"         ,
	"cat"          => "Cat"           ,
	"dog"          => "Dog"           ,
	"microbat"     => "Microbat"      ,
	"megabat"      => "Megabat"       ,
	"hedgehog"     => "Hedgehog"      ,
	"shrew"        => "Shrew"         ,
	"elephant"     => "Elephant"      ,
	"rock_hyrax"   => "Rock_hyrax"    ,
	"tenrec"       => "Tenrec"        ,
	"armadillo"    => "Armadillo"     ,
	"sloth"        => "Sloth"
	);
}
else{
    %phylogenyNames =(
	"dmel" => "dmel" ,
	"dsim" => "dsim" ,
	"dsec" => "dsec" ,
	"dyak" => "dyak" ,
	"dere" => "dere" ,
	"dana" => "dana" ,
	"dpse" => "dpse" ,
	"dper" => "dper" ,
	"dwil" => "dwil" ,
	"dvir" => "dvir" ,
	"dmoj" => "dmoj" ,
	"dgri" => "dgri"
	);
}

#make fasta_aln
my $cmd_fasta_aln = "ls $outAlignments | xargs -I X sh -c \'t_coffee -other_pg seq_reformat -in $outAlignments/X -output fasta_aln > ${fasta_aln}/X\'";
(system "$cmd_fasta_aln") == 0 or die "Error[csf.pl]! cannot execute\n$cmd_fasta_aln\n$!\n";



if (-f "${pipelineDirName}/experiments/${experimentName}/csfOut.txt" ) {system "rm ${pipelineDirName}/experiments/${experimentName}/csfOut.txt";}


#convertNames
opendir (FA,"${fasta_aln}") or die "Error[csf.pl]! Cannot open the ${fasta_aln} dir \n$!\n";
my @fasta_aln_files = readdir(FA);
closedir(FA);
my @msaSpeciesNames;
foreach my $f (@fasta_aln_files){
  next if (($f eq '.') or ($f eq '..'));
  my $fasta_aln_file = "${fasta_aln}/$f";
  open (F,"<$fasta_aln_file") or die "Error[csf.pl]! Cannot open the $fasta_aln_file file \n$!\n";
  my $c   = 0;
  my $spy = 0;
  foreach my $line (<F>){
    $c++;
    chomp $line;

    #reference species is always the first sequence
    if ($c == 1) {
      system "sed -i \'${c}c\\>$phylogenyNames{$referencePhyloCsfName}\' $fasta_aln_file";
      $spy = 2;
      push (@msaSpeciesNames , $phylogenyNames{$referencePhyloCsfName});
      next;
    }

    if (($spy == 1) and ($line!~/^>(\S+)/)){system "sed -i \'${c}c\\\' $fasta_aln_file" ; $c-- ; next ;}

    if ($line=~/^>(\S+)/){
      $spy = 1;
      my $rowHeader = lc($1);
      foreach my $name (keys %phylogenyNames){
	next if ($name eq lc ($referencePhyloCsfName)); #skipp the reference species duplicate
	if ($rowHeader =~/$name/){
	  system "sed -i \'${c}c\\>$phylogenyNames{$name}\' $fasta_aln_file";
	  $spy = 2;
	  push (@msaSpeciesNames , $phylogenyNames{$name});
	}
      }
      if ($spy == 1){
	system "sed -i \'${c}c\\\' $fasta_aln_file";
	$c--;
      }
    }

  }
  close F;


  #skipp
  if (scalar @msaSpeciesNames < 2 ){system "echo \"$fasta_aln_file NA\">> ${pipelineDirName}/experiments/${experimentName}/csfOut.txt"; next;}

  #RUN pyloCSF
  my $cmd_phyloCSF = "$phylocsfName $phyloCSFparametersName $fasta_aln_file --removeRefGaps --species=" ;
  my $i = 0;
  foreach my $s (@msaSpeciesNames){
    $cmd_phyloCSF .= $s;
    $i++;
    $cmd_phyloCSF .= "," if ($i <= scalar @msaSpeciesNames);
  }
  $cmd_phyloCSF .= " >> ${pipelineDirName}/experiments/${experimentName}/csfOut.txt 2>> ${pipelineDirName}/experiments/${experimentName}/csfOut.txt";
  (system "$cmd_phyloCSF") == 0 or die "Error[csf.pl]! cannot run $cmd_phyloCSF \n$!\n";
  @msaSpeciesNames = ();
}








###
#FUNCTIONS
sub options {
  my ($tcoffeeName  , $phylocsfName , $experimentName , $pipelineDirName , $phyloCSFparametersName , $referencePhyloCsfName);
  my $spyReferencePhyloCsf = 1;
  my $spyPhylocsfParameter = 1;
  my $spyPhylocsf          = 1;
  my $spyTcoffee           = 1;
  my $spyExperiment        = 1;
  my $spyPipelineDir       = 1;
  foreach my $field (0..$#ARGV){
    if ($ARGV[$field] eq '-phyloCSFreference'){
	$referencePhyloCsfName = $ARGV[1+$field];
	$spyReferencePhyloCsf = 2;
	next;
    }
    if ($spyReferencePhyloCsf == 2){
	$spyReferencePhyloCsf = 3;
	next;
    }
    if ($ARGV[$field] eq '-phyloCSFparameters'){
	$phyloCSFparametersName = $ARGV[1+$field];
	$spyPhylocsfParameter = 2;
	next;
    }
    if ($spyPhylocsfParameter == 2){
	$spyPhylocsfParameter = 3;
	next;
    }
    if ($ARGV[$field] eq '-phyloCSF'){
	$phylocsfName = $ARGV[1+$field];
	$spyPhylocsf = 2;
	next;
    }
    if ($spyPhylocsf == 2){
	$spyPhylocsf = 3;
	next;
    }
    if ($ARGV[$field] eq '-t_coffee'){
	$tcoffeeName = $ARGV[1+$field];
	$spyTcoffee = 2;
	next;
    }
    if ($spyTcoffee == 2){
	$spyTcoffee = 3;
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
  die "Error[csf.pl]! You must provide the -experiment parameter.\n"     if ($spyExperiment  != 3);
  die "Error[csf.pl]! You must provide the -pipeline_dir parameter.\n"   if ($spyPipelineDir != 3);
  if ((defined $phyloCSFparametersName) && ($phyloCSFparametersName ne '29mammals') && ($phyloCSFparametersName ne '12flies')){
    die "Error[csf.pl]! You must provide the -phyloCSFparameters field with \'29mammals\' or \'12flies\'\n";
  }
  $phyloCSFparametersName = '29mammals'  if (! defined $phyloCSFparametersName);
  $tcoffeeName            = "t_coffee"   if (! defined $tcoffeeName);
  $phylocsfName           = "PhyloCSF"   if (! defined $phylocsfName);
  $referencePhyloCsfName  = "human"      if (! defined $referencePhyloCsfName);

  return ($tcoffeeName  , $phylocsfName , $experimentName , $pipelineDirName , $phyloCSFparametersName , $referencePhyloCsfName);
}
sub acceptedVariableSpace {
  my %space = ('-t_coffee' => 1 , '-phyloCSF' => 1  , '-experiment' => 1  , '-pipeline_dir' => 1 , '-phyloCSFparameters' => 1 , '-phyloCSFreference' => 1);
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[csf.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}


sub help_message {
my $helpMessage = "\nNAME
csf.pl - Convert the clustal MSA into fasta_aln format. Then it runs PhyloCSF to estimate the codon substitution frequency (CSF)  \n
SYNOPSIS
csf.pl -experiment -pipeline_dir [-t_coffee -phyloCSF -phyloCSFparameters -phyloCSFreference] \n
DESCRIPTION
   * csf.pl takes as inputjust the path to tcoffee, phyloCSF, and the experiment and the pipeline dir name to find the alignments
   * it fills the \"fasta_aln\" in the experiment folder and it returns a file with the CSF scores associated to each MSA.
OPTIONS
   * By default csf.pl assumes that tcoffee program name on your computer is \'t_coffee\'.
     If this is not the case you must provide the field -t_coffee with the appropriate path.
   * By default csf.pl assumes that pyloCSF program name on your computer is \'pyloCSF\'.
     If this is not the case you must provide the field -pyloCSF with the appropriate path.
   * -phyloCSFparameters [29mammals|12flies]   default 29mammals
   * -phyloCSFreference  [species Name]        default human



TROUBLESHOOTING
   * pyloCSF so far can be used just on two different phylogenies, 29 mammals and 12 flies. If your data do not belong to one of these phylogenies (for at least two species) csf.pl cannot be executed.

";
print "$helpMessage\n\n\n";
exit;
}



