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
my ($queryName , $experimentName , $pipelineDirName) = options();
my $dataDir  = "${pipelineDirName}/experiments/${experimentName}/outSim";

my @allSpecies = `ls ${pipelineDirName}/allGenomeInfo/`; if ($?){die "Error[sim2matrix.pl]! Cannot do:\nls ${pipelineDirName}/allGenomeInfo/\n$!\n";}
print "species,";
my $i = 0;
foreach my $species (@allSpecies){
  next if (($species eq '.') or ($species eq '..'));
  $i++;
  chomp $species ;
  print "$species";
  print "," if ($i < scalar(@allSpecies));
}
print "\n";


my @notFound;
open (Q,"<$queryName") or die "Error[sim2matrix.pl]!";
foreach my $line (<Q>){
    chomp $line;
    if ($line =~ /^>(.+)/){
	my $id = $1;
	push (@notFound , $id) if (! -e "${dataDir}/$id");
    }
}
close Q;

opendir (S,"$dataDir") or die "Error[sim2matrix.pl]! Cannot read $dataDir \n$!\n";
my @simFiles = readdir(S);
closedir S;

foreach my $tx_id (@simFiles){
    next if (($tx_id eq '.') or ($tx_id eq '..'));
    my %info;
    my $tx_sim_file = $dataDir . "/$tx_id";

    open (SIM,"<$tx_sim_file") or die "Error[sim2matrix.pl]! Cannot open the sim file $tx_sim_file\n";
    foreach my $line (<SIM>){
      chomp $line;
      #initialize the names
      if ($line=~/^# SEQ_INDEX (\S+) (\S+)/){
	my $id    = $1;
	my $index = $2;
	if ($id !~/_hit/){
	  die "Error[sim2matrix.pl]!the id is $id and the $index is $index...error\n" if ($index != 0);
	  next;
	}
	else{
	  if ($id =~/hit\d+_(\S+)$/){
	    $id = $1;
	    $info{$tx_id}{$id} = 0;
	  }
	}
      }
      if ($line=~/BOT\s+\d+\s+\d+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)/){
	my $currentID1 = $1;
	my $currentID2 = $2;
	my $score      = $3;

	#take species name
	if ($currentID2 =~/hit\d+_(\S+)$/){
	  $currentID2 = $1;
	}

	#find the comparison between the original query and each homolog
	if ($currentID1 eq $tx_id){
	  #if there are multiple homologs, just take the one with highest percent id
	  if(! defined $info{$tx_id}{$currentID2}){
	    $info{$tx_id}{$currentID2} = $score;
	  }
	  else{
	    $info{$tx_id}{$currentID2} = max2($info{$tx_id}{$currentID2} , $score);
	  }
	}

      }
    }
    close SIM;



    my $i = 0;
    print "$tx_id,";
    foreach my $species (@allSpecies){
	$i++;
	if (defined $info{$tx_id}{$species}){
	    print $info{$tx_id}{$species} ;
	}
	else{
	    print "0";
	}
	print "," if ($i < scalar(@allSpecies));
    }
    print "\n";

}



if (@notFound){
    foreach my $tx_id (@notFound){
	my $i = 0;
	print "$tx_id,";
	foreach my $species (@allSpecies){
	    $i++;
	    print "0";
	    print "," if ($i < scalar(@allSpecies));
	}
	print "\n";
    }
}





#FUNZIONI
sub options {
  my ($queryName , $experimentName , $pipelineDirName);
  my $spyExperiment                = 1;
  my $spyPipelineDir               = 1;
  my $spyQuery                     = 1;


  foreach my $field (0..$#ARGV){
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
  die "Error[sim2matrix.pl]! Please provide the -query parameter\n"   if (! defined $queryName);
  die "Error[sim2matrix.pl]! Please provide the -experiment parameter\n"   if (! defined $experimentName);
  die "Error[sim2matrix.pl]! Please provide the -pipeline_dir parameter\n" if (! defined $pipelineDirName);
  return ($queryName , $experimentName , $pipelineDirName);
}

sub acceptedVariableSpace {
  my %space = ('-pipeline_dir' => 1 , '-experiment' => 1 , '-query' => 1);
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[sim2matrix.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}

sub help_message {
my $helpMessage = "\nNAME
sim2matrix.pl - It generates the simMatrix.csv out of the outSim directory\n
SYNOPSIS
sim2matrix.pl -pipeline_dir -experiment -query 

DESCRIPTION
   * sim2matrix.pl creates a comma separated values matrix (.csv) out of considering the similarity scores among the detected homologs
   * Transcripts that were included in the query file but are not found in the outSim directory will be anyway reported and represented with zeros
   * These are in fact transcripts that failed either the blastSearch either the exonerate remapping steps.

OPTIONS
   * This script does not accept any option

";
print "$helpMessage\n\n\n";
exit;
}

sub max2 {
  my ($a, $b) = @_;

  return $b unless (defined $a);
  return $a unless (defined $b);
  return (($a > $b)? $a: $b);
}
