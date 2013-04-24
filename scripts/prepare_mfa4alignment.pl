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
my ($query  ,  $experimentName , $pipelineDirName) = options();
my $exonerate_out = "${pipelineDirName}/experiments/${experimentName}/EXONERATE_OUT/";
my $prepareOutDir = "${pipelineDirName}/experiments/${experimentName}/multifasta4EachTx/";


#Read the multifasta
open (M,"<$query") or die "Error[prepare_mfa4alignment.pl]!cannot read $query $!\n";
my @transcripts = <M>;
close M;
my %allQueries =  loadMultifastaIntoHash (@transcripts);

#create a file for each query
print "#Creating a homolog multifasta FASTA file for each query...\n";
initiateMFA();



#find and concatenate the orthologs
opendir ( DIR, $exonerate_out ) || die "Error[prepare_mfa4alignment.pl]! cannot open the directory $exonerate_out\n";
my @allSpeciesOrthologhFiles = readdir(DIR);
closedir(DIR);

foreach my $file (@allSpeciesOrthologhFiles){
  next if (($file eq '.') or  ($file eq '..') or ($file=~/~/));
  if ($file=~/\.gtf$/){next;}
  my $fileName    = "$exonerate_out" . "/$file";
  my $speciesName;
  if ($file =~/^(.+)\.fa$/){
    $speciesName = $1;
  }
  else {die "Error[prepare_mfa4alignment.pl]! Strange exonerate output file name:\n$fileName\n";}

  foreach my $query (keys %allQueries){
    my $queryName = $query;
    $queryName =~s/^>//;

    my @fastaFetcherOut = fastaFetch ($fileName , $queryName);

    next if (! defined $fastaFetcherOut[0]);

    foreach my $l (0..$#fastaFetcherOut){
      if ($fastaFetcherOut[$l] =~/^>/){
	chomp $fastaFetcherOut[$l];
	$fastaFetcherOut[$l] .= "_$speciesName";
	$fastaFetcherOut[$l] .= "\n";
      }
    }
    open (O,">>$prepareOutDir/${queryName}.mfa");
    foreach my $l (@fastaFetcherOut){print O "$l";}
    close O;
  }


}


#FUNCTIONS
sub loadMultifastaIntoHash{
  my (@multifasta) = @_;
  my (%allTheSequences , $sequence , $header);
  my $spy = 1;
  foreach my $line (@multifasta){
    chomp $line;
    next if ($line=~/^\s*$/);
    if ($line=~/(>\S+)/){
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

sub help_message {
my $helpMessage = "\nNAME
prepare_mfa4alignment.pl - Starting from the exonerateRemapping.pl output, it returns a .mfa file for each query transcript

SYNOPSIS
prepare_mfa4alignment.pl -query file -experiment -pipeline_dir

DESCRIPTION
   * prepare_mfa4alignment.pl takes the output returned by exonerateRemapping.pl and the original multi-FASTA query file
   * It returns in \"multifasta4EachTx\" directory the .mfa files embedding the homologs detected for each query.

";
print "$helpMessage\n\n\n";
exit;
}


sub options {
  my ($query , $experimentName , $pipelineDirName);
  my $spyQuery                = 1;
  my $spyExperiment           = 1;
  my $spyPipelineDir          = 1;

  foreach my $field (0..$#ARGV){
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
   if ($ARGV[$field] eq '-query'){
      $query = $ARGV[1+$field];
      $spyQuery = 2;
      next;
    }
    if ($spyQuery == 2){
      $spyQuery = 3;
      next;
    }

  }
  die "Error[prepare_mfa4alignment.pl]! You must provide the -query parameter. Please indicate the multi-FASTA query file\n" if ($spyQuery != 3);
  die "Error[prepare_mfa4alignment.pl]! Please provide the -experiment parameter\n"                                          if (! defined $experimentName);
  die "Error[prepare_mfa4alignment.pl]! Please provide the -pipeline_dir parameter\n"                                        if (! defined $pipelineDirName);
  return ($query , $experimentName , $pipelineDirName);
}

sub acceptedVariableSpace {
  my %space = ('-experiment' => 1 , '-pipeline_dir' => 1   , '-query' => 1);
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[prepare_mfa4alignment.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}

sub fastaFetch {
    my ($fileName , $nameToExtract) = @_;
    #putting all the FASTA sequences in an hash of arrays
    my $spy = 0;
    my ($name , %allFASTA ,  @headers );
    open (FILE1,"<$fileName") or die ("Error[prepare_mfa4alignment]!. Cannot open the multi fasta file:\n$fileName\n$!\n");
    foreach my $line (<FILE1>){
	chomp ($line);
	if ($line=~/^>(.+)$/){
	    $name = $1;
	    $spy  = 1;
	    push (@{$allFASTA{$name}},"$line\n");
	    push (@headers , "$name");
	    next;
	}
	if (($spy == 1) && ($line=~/\w+/)){
	    push (@{$allFASTA{$name}},"$line\n");
	    next;
	}
	#return space
	if (($spy == 1) && ($line=~/^\s*$/)){
	    $spy = 0;
	    next;
	}
    }
    close (FILE1);

    #extraction
    my ( %matches , @out );
    chomp $nameToExtract;
    foreach my $header (@headers){
	if ($header=~/$nameToExtract/){
	    $matches{$header} = 1;
	    foreach my $line (@{$allFASTA{$header}}) {
		push (@out , "$line");
	    }
	}
      }
    return @out;
}


sub initiateMFA {
  foreach my $query (keys %allQueries){
    my $queryName = $query;
    $queryName =~s/^>//;
    my $sequence = string2FASTA($allQueries{$query});
    open (Q,">$prepareOutDir/${queryName}.mfa") or die "cannot open $queryName $!\n";
    print Q "$query\n$sequence\n";
    close Q;
  }
}
