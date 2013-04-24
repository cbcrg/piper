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
my ($queryName , $referenceName , $experimentName , $pipelineDirName) = options();
my $dataDir  = "${pipelineDirName}/experiments/${experimentName}/EXONERATE_OUT";
my $outDir  = "${pipelineDirName}/experiments/${experimentName}/results";

system "mkdir -p $outDir";
open (O1,">${outDir}/familyMembers.csv") or die "Error[family2matrix.pl]! Cannot create ${outDir}/familyMembers.csv $!\n";
open (O2,">${outDir}/familyExpansionContraction.csv") or die "Error[family2matrix.pl]! Cannot create ${outDir}/familyExpansionContraction.csv $!\n";

#take species name
my @allSpecies = `ls ${pipelineDirName}/allGenomeInfo/`; if ($?){die "Error[family2matrix.pl]! Cannot do:\nls ${pipelineDirName}/allGenomeInfo/\n$!\n";}



#take the number of homologs in each species
my (%info , %info2 , @allIds);
open (Q,"<$queryName") or die "Error[family2matrix.pl]!";
while (my $line =<Q>){
  chomp $line;
  if ($line =~ /^>(.+)/){
    my $id = $1;
    push(@allIds , $id);
    foreach my $species (@allSpecies){
      chomp $species ;
      next if (($species eq '.') or ($species eq '..'));
      my $number = `grep -cP \">${id}_hit*\" ${dataDir}/${species}.fa | tr -d \"\n\"`;
      $info{$id}{$species} = $number;
    }
  }
}
close Q;


#take the difference in homologs with respect to the reference species
foreach my $id (@allIds){
  my $referenceNumber =  $info{$id}{$referenceName};
  foreach my $species (@allSpecies){
    next if (($species eq '.') or ($species eq '..'));
    $info2{$id}{$species} = $info{$id}{$species} - $referenceNumber;
  }
}


#printOut1
my $i = 0;
print O1 "species,";
foreach my $species (@allSpecies){
  next if (($species eq '.') or ($species eq '..'));
  $i++;
  chomp $species ;
  print O1 "$species";
  print O1"," if ($i < scalar(@allSpecies));
}
print O1 "\n";

foreach my $id (@allIds){
  my $i = 0;
  print O1 "$id,";
  foreach my $species (@allSpecies){
    chomp $species ;
    next if (($species eq '.') or ($species eq '..'));
    $i++;
    print O1 $info{$id}{$species};
    print O1 "," if ($i < scalar(@allSpecies));
  }
  print O1 "\n";
}

#printOut2
$i = 0;
print O2 "species,";
foreach my $species (@allSpecies){
  chomp $species ;
  next if (($species eq '.') or ($species eq '..'));
  $i++;
  chomp $species ;
  print O2 "$species";
  print O2 "," if ($i < scalar(@allSpecies));
}
print O2 "\n";

foreach my $id (@allIds){
  my $i = 0;
  print O2 "$id,";
  foreach my $species (@allSpecies){
    next if (($species eq '.') or ($species eq '..'));
    $i++;
    print O2 $info2{$id}{$species};
    print O2 "," if ($i < scalar(@allSpecies));
  }
  print O2 "\n";
}

close O1;
close O2;








#FUNZIONI
sub options {
  my ($queryName , $referenceName , $experimentName , $pipelineDirName);
  my $spyExperiment                = 1;
  my $spyPipelineDir               = 1;
  my $spyQuery                     = 1;
  my $spyReferenceName             = 1;


  foreach my $field (0..$#ARGV){
    if ($ARGV[$field] eq '-referenceName'){
	$referenceName = $ARGV[1+$field];
	$spyReferenceName = 2;
	next;
    }
    if ($spyReferenceName == 2){
	$spyReferenceName = 3;
	next;
    }
    if ($ARGV[$field] eq '-query'){
	$queryName = $ARGV[1+$field];
	$spyQuery = 2;
	next;
    }
    if ($spyQuery == 2){
	$spyQuery = 3;
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
  die "Error[family2matrix.pl]! Please provide the -query parameter\n"         if (! defined $queryName);
  die "Error[family2matrix.pl]! Please provide the -experiment parameter\n"    if (! defined $experimentName);
  die "Error[family2matrix.pl]! Please provide the -pipeline_dir parameter\n"  if (! defined $pipelineDirName);
  die "Error[family2matrix.pl]! Please provide the -referenceName parameter\n" if (! defined $referenceName);

  return ($queryName , $referenceName ,$experimentName , $pipelineDirName);
}

sub acceptedVariableSpace {
  my %space = ('-pipeline_dir' => 1 , '-referenceName' => 1 , '-experiment' => 1 , '-query' => 1);
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[family2matrix.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}

sub help_message {
my $helpMessage = "\nNAME
family2matrix.pl - It generates the familyMatrix.csv out of the EXONERATE_OUT directory\n
SYNOPSIS
family2matrix.pl -pipeline_dir -experiment -query -referenceName

DESCRIPTION
   * family2matrix.pl creates a comma separated values matrix (.csv) out of considering family expansion and family restriction events
   * The referenceName is the name of the reference species with respect to whom the events of family expansion and family restriction are considered (e.g. human)
   * The query is the fasta file that was used by pipeR. Tipically the queries are transcripts coming from the reference species.

OPTIONS
   * This script does not accept any option

";
print "$helpMessage\n\n\n";
exit;
}
