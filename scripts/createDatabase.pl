#!/usr/bin/env perl
use strict;
use warnings;
#Author: Giovanni Bussotti

#HELP
my $help  = 0;
my $doNotSplitReferenceGenome = 0;
foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '-help') or ($ARGV[$field] eq '--help'));
}
if ($help > 0){
 help_message();
}
#TAKE OPTIONS
acceptedVariableSpace();
my ($experimentName , $referenceGenomeName , $splitGenome, $xdFormatName  , $inputGenomes , $strategyName , $pipelineDirName) = options();

#READ INPUT GENOMES
my @allMFAfiles;
if (-f $inputGenomes){
  open (F,"<$inputGenomes") || die "Error[createDatabase.pl]!cannot read the $inputGenomes file $!\n";
  @allMFAfiles = <F>;
  close F;
}
else {die "Error[createDatabase.pl]! $inputGenomes does not exist $!\n";}


#CREATE OUT and STDERR DIR
my $outDirName = $pipelineDirName . "/allGenomeInfo" ;


print "#formatting ..\n";
foreach my $line (@allMFAfiles){
    chomp $line;
    next if ($line=~/^\s*$/);
    my ($genomeName , $inputFileName);
    if ($line=~/^(\S+)\s+(\S+)/){
      $inputFileName = $1;
      $genomeName    = $2;
    }
    else{ die "Error[createDatabase.pl]! in $inputGenomes the line:\n$line\nmust have this syntax:\nfileName outputGenomeName\nPlease fix it and re-run me!\n" ;}
    die "Error[createDatabase.pl]! $inputFileName included in $inputGenomes does not exist. please specify the FULL file path\n" unless (-f $inputFileName);
    my $outGenomeDirName = $outDirName . "/$genomeName";

    #SKIPP ALREADY EXISTING GENOMES
    if ((-d $outGenomeDirName) && (`ls -A $outGenomeDirName`)) {next;}
    (system "mkdir -p $outGenomeDirName") == 0 or die "Error[createDatabase.pl]! Cannot create $outGenomeDirName genome folder \n$!\n";

    #MAKE CHROMOSOME DIR
    (system "mkdir ${outGenomeDirName}/chr") == 0 or die "Error[createDatabase.pl]! Cannot create ${outGenomeDirName}/chr genome folder \n$!\n";
    my $sanityCheck = multifastaSplitter ($inputFileName , "${outGenomeDirName}/chr");
    if ($sanityCheck == 1){
      system "rm -rf ${outGenomeDirName}";
      exit 1;
    }
    #MAKE DATABASE
    makeDatabase($outGenomeDirName , $inputFileName , $genomeName);



    print "#$genomeName\n";
}


#PREPARE REFERENCE DATABASE
if ($referenceGenomeName ne 'none'){
  unless (-d "$pipelineDirName/experiments/$experimentName/pre_processing/inputInfo/allGenomeInfo/"){
    my $inputInfoDir = "$pipelineDirName/experiments/$experimentName/pre_processing/inputInfo/";
    (system "mkdir $inputInfoDir")                                == 0 or die "Error[orthologFinder]! cannot create the inputInfo Folder\n$!\n";
    (system "mkdir  ${inputInfoDir}/allGenomeInfo")               == 0 or die "Error[orthologFinder]! cannot create the allGenomeInfo Folder\n$!\n";
    (system "mkdir  ${inputInfoDir}/allGenomeInfo/reference")     == 0 or die "Error[orthologFinder]! cannot create the reference Folder\n$!\n";
    (system "mkdir  ${inputInfoDir}/allGenomeInfo/reference/chr") == 0 or die "Error[createDatabase.pl]! Cannot create ${inputInfoDir}/allGenomeInfo/chr reference genome folder \n$!\n";
    $doNotSplitReferenceGenome = 1;
    my $sanityCheck = multifastaSplitter ($referenceGenomeName , "${inputInfoDir}/allGenomeInfo/reference/chr");
    makeDatabase( "${inputInfoDir}/allGenomeInfo/reference" , $referenceGenomeName , 'referenceGenome');
    print "#formatting reference\n";
  }
}



###
#FUNCTIONS
sub help_message {
my $helpMessage = "\nNAME
createDatabase.pl - Input a file specifying a list of multi-FASTA (absolute path) to have the RNAmappingPipeline database in the \"allGenomeInfo\" pipeline directory\n
SYNOPSIS
createDatabase.pl -genomes -pipeline_dir [-strategy -xdformat]\n
DESCRIPTION
   * CreateDatabase.pl takes as input a file listing the genomes with this syntax:  genomeFileName(absolute path) genomeName
   * For \"genomeFileName\", this script will return a \"genomeName\" folder including the xdformat database and a \"chr\" folder encompassing the chromosome FASTA files.
OPTIONS
   * By default createDatabase.pl assumes that xdformat program name on your computer is just \'xdformat\'. If this is not the case you must provide the field -pg with the xdformat path.
   * By default createDatabase.pl assumes that the blast flavor that will be used is wublastn. You can change this by editing the -strategy parameter
TROUBLESHOOTING
   * The \"allGenomeInfo\" is unique across all the experiment in the pipeline directory.
     Anytime the user start a new experiment createDatabase.pl checks that all the genomes specified in the new genomes file (-genomes) already exist in \"allGenomeInfo\", if not it creates them.
     All the genomes that already existed will be kept.
   * Both AB-Blast and WU-Blast packages are provided with their own xdformat. If needed it is therefore important that the user properly specifies, with -pg option, the proper one. The same for blastr.
     Mind to use the xdformat included in the blast package you want to use
";
print "$helpMessage\n\n\n";
exit;
}

sub makeDatabase {
    my ($outGenomeDirName , $inputFileName , $genomeName)= @_;
    my $db_folder; 
    if (($strategyName eq 'wublastn') or ($strategyName eq 'wublastn_opt')){$db_folder = "$outGenomeDirName".'/wublastn_db';}
    if (($strategyName eq 'abblastn') or ($strategyName eq 'abblastn_opt')){$db_folder = "$outGenomeDirName".'/abblastn_db';}
    if ($strategyName eq 'wublastr'){$db_folder = "$outGenomeDirName".'/wublastr_db';}
    if ($strategyName eq 'abblastr'){$db_folder = "$outGenomeDirName".'/abblastr_db';}
    (system "mkdir $db_folder" ) == 0 or die "Error[createDatabase.pl]! Cannot create the $db_folder folder $!\n";

    if (($splitGenome eq 'yes') && ($doNotSplitReferenceGenome == 0)) {
      my $splitGenomeCmd = "for X in `ls ${outGenomeDirName}/chr/`; do $xdFormatName -n -o ${db_folder}/\$X - < ${outGenomeDirName}/chr/\$X ; done";
      (system "$splitGenomeCmd") == 0 or die "Error[createDatabase.pl]! Cannot execute the splitGenome command\n$splitGenomeCmd\n$!\n";
      return;
    }
    my @names = `ls --color=never ${outGenomeDirName}/chr/`; if ($?) {die "Error[createDatabase.pl]! cannot take the chr file names to format the database $!\n"};
    my $format_cmd;
    if (scalar @names > 250){
      print STDERR "#Genome $genomeName has more than 250 scaffolds. The formatting will take a bit of time..\n";
      my $tmp_genome = fileNameGenerator("$pipelineDirName/experiments/$experimentName/$genomeName");
      foreach my $c (@names){
	chomp $c;
	(system "cat ${outGenomeDirName}/chr/$c >> $tmp_genome") == 0 or die "Error[createDatabase.pl]! Cannot create $tmp_genome\n$!\n";
      }
      $format_cmd = "$xdFormatName -n -o ${db_folder}/db $tmp_genome";
      (system "$format_cmd") == 0 or die "Error[createDatabase.pl]! Cannot run $format_cmd\n$!\n";;
      system "rm $tmp_genome";
      return;
    }
    $format_cmd = "$xdFormatName -n -o ${db_folder}/db - <";
    foreach my $name (@names){ chomp $name; $format_cmd .= " ${outGenomeDirName}/chr/$name";}
    (system "$format_cmd") == 0 or die "Error[createDatabase.pl]! Cannot run $format_cmd\n$!\n";;
  }

sub fileNameGenerator{
  my ($nameRoot) = @_;
  my $tmp_name_counter = 0;
  my $tmp_name;
  while (!$tmp_name || -f $tmp_name) {
    $tmp_name_counter++;
    $tmp_name = "${nameRoot}_$$".".${tmp_name_counter}";
  }
  return $tmp_name;
}

sub options {
  my ($experimentName , $referenceGenomeName , $splitGenome , $xdformatName  , $genomesName , $strategyName , $pipelineDirName);
  my $spyXdformat        = 1;
  my $spyGenomes         = 1;
  my $spyStrategy        = 1;
  my $spyPipelineDir     = 1;
  my $spySplitGenome     = 1;
  my $spyReferenceGenome = 1;
  my $spyExperiment      = 1;

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
    if ($ARGV[$field] eq '-reference_genome'){
      $referenceGenomeName = $ARGV[1+$field];
      $spyReferenceGenome = 2;
      next;
    }
    if ($spyReferenceGenome == 2){
      $spyReferenceGenome = 3;
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
    if ($ARGV[$field] eq '-splitGenome'){
      $splitGenome = $ARGV[1+$field];
      $spySplitGenome = 2;
      next;
    }
    if ($spySplitGenome == 2){
      $spySplitGenome = 3;
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
    if ($ARGV[$field] eq '-genomes'){
      $genomesName = $ARGV[1+$field];
      $spyGenomes = 2;
      next;
    }
    if ($spyGenomes == 2){
      $spyGenomes = 3;
      next;
    }
    if ($ARGV[$field] eq '-strategy'){
	$strategyName = $ARGV[1+$field];
	$spyStrategy = 2;
	next;
    }
    if ($spyStrategy == 2){
	$spyStrategy = 3;
	next;
    }
  }

  die "Error[createDatabase.pl]! You must provide the -genomes parameter. Please indicate either a folder containing the input genomes, either a file listing the pointers to the genome files\n" if ($spyGenomes != 3);
  die "Error[createDatabase.pl]! You must provide the -pipeline_dir parameter.\n" if ($spyPipelineDir != 3);

  $xdformatName     = 'xdformat' if (! defined $xdformatName);
  if ((defined $strategyName) && ($strategyName ne "wublastn") && ($strategyName ne "abblastn") && ($strategyName ne "wublastr") && ($strategyName ne "abblastr") && ($strategyName ne "abblastn_opt") && ($strategyName ne "wublastn_opt")){
      die "Error[createDatabase.pl]! The supported blast flavors are wublastn abblastn blastr abblastn_opt wublastn_opt\n";
  }
  die "Error[createDatabase.pl]! -splitGenome parameter can be either \'yes\' or \'no\'\n" if ((defined $splitGenome) && ($splitGenome ne 'yes') && ($splitGenome ne 'no'));
  $strategyName        = 'wublastn'  if (! defined $strategyName);
  $splitGenome         = 'no'        if (! defined $splitGenome);
  $referenceGenomeName = 'none'      if (! defined $referenceGenomeName);
  die "Error[createDatabase.pl]! If you specify the reference_genome then you need to tell also what is the experiment name\n"   if ((! defined $experimentName) && ($referenceGenomeName ne 'none'));
  return ($experimentName , $referenceGenomeName , $splitGenome , $xdformatName , $genomesName , $strategyName , $pipelineDirName);
}

sub acceptedVariableSpace {
  my %space = ('-experiment' => 1  , '-reference_genome' => 1 , '-splitGenome' => 1 , '-pipeline_dir' => 1  , '-xdformat' => 1  , '-genomes' => 1  , '-strategy' => 1);
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[createDatabase.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
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

sub multifastaSplitter {
  my ($multiFastaFile , $outDir) = @_;
  open (I,"<$multiFastaFile") or die "Error[createDatabase.pl]! cannot open the input multifasta file $!\n";

  #SPLITTING
  my $nameSpy = 0;
  my $jolly   = 0;
  my ($name , $fileName);
  while (my $line = <I>){
    chomp $line;
    next if ($line=~/^\s*$/);

    if ($line=~/^>(\S+)/){
      $name     = $1;
      $name =~s/[\|| |\>|\<|\/|\'|\:|\\]/_/g;
      $fileName = $name;

      if ($nameSpy == 1){
	$nameSpy = 0;
	$jolly   = 0;
	close F;
      }
      open (F,">$outDir"."/$fileName") or die "Error[createDatabase.pl]!cannot create the $outDir"."/$fileName FASTA file $!\n";
      #$name  =~ s/[\|| |\>|\<]/_/g;
      print F ">$name\n";
      $jolly = 1;
      next;
    }

    if (($line!~/>/) and ($jolly == 1)){
      if (length($line) > 100){
	print F string2FASTA($line) . "\n";
      }
      else{
	print F "$line\n";
      }
      $nameSpy = 1;
      next;
    }

  }
  close F;
  close I;
  return 0;
}







###############OLD
###This function is perfectly working. However, to make the pipeline compatible with the Server we decided not to have a strict sanity check of the taget genome headers, but just take everything before the first space.
# sub multifastaSplitter {
#   my ($multiFastaFile , $outDir) = @_;
#   open (I,"<$multiFastaFile") or die "Error[createDatabase.pl]! cannot open the input multifasta file $!\n";
#   #SPLITTING
#   my $nameSpy = 0;
#   my $jolly   = 0;
#   my ($name , $fileName);
#   foreach my $line (<I>){
#     chomp $line;
#     next if ($line=~/^\s*$/);
#     if ($line=~/^>(.+)$/){
#       $name     = $1;
#       if ($name=~/[\|| |\>|\<|\/|\'|\:]/){
# 	print "Error[createDatabase.pl]! the chromosome (or contig) header $name of $multiFastaFile contains one of the following forbidden symbol: \'\:\' , \'\'\' ,  \'\|\'  ,  \' \'  ,  \'\>\'  ,  \'\<\'  ,  \'\/\'\nPlease replace it with something else like \'_\', and re-run the script\n";
# 	close F if ($nameSpy == 1);
# 	close I;
# 	return 1;
#       }
#       $fileName = $name;
#       #$fileName =~ s/\//-/g;
#       if ($nameSpy == 1){
# 	$nameSpy = 0;
# 	$jolly   = 0;
# 	close F;
#       }
#       open (F,">$outDir"."/$fileName") or die "Error[createDatabase.pl]!cannot create the $outDir"."/$fileName FASTA file $!\n";
#       #$name  =~ s/[\|| |\>|\<]/_/g;
#       print F ">$name\n";
#       $jolly = 1;
#       next;
#     }
#     if (($line!~/>/) and ($jolly == 1)){
#       if (length($line) > 100){
# 	print F string2FASTA($line) . "\n";
#       }
#       else{
# 	print F "$line\n";
#       }
#       $nameSpy = 1;
#       next;
#     }
#   }
#   close F;
#   close I;
#   return 0;
# }
























