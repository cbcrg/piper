#!/usr/bin/env perl
use strict;
use warnings;
#Author: Giovanni Bussotti

#HELP
my $help  = 0;
my ($cmd , %infoGtf);
foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '-help') or ($ARGV[$field] eq '--help'));
}
if ($help > 0){
 help_message();
}
#TAKE OPTIONS
acceptedVariableSpace();
my ($splitGenome , $unstrandedGTFname, $link2referenceGenomeName , $xdformatName , $gtfName , $fastaName , $referenceGenomeName , $chr_subseqName , $clusterName  , $strategyName , $blastName , $exonerate_success_mode , $exonerate_lines_mode  , $experimentName , $pipelineDirName ) = options();
my $clusterConfigName          = "${pipelineDirName}/experiments/${experimentName}/CONFIG/clusterConfig";
my $exonerateExtensionFileName = "${pipelineDirName}/experiments/${experimentName}/CONFIG/exonerateExtensionFile";
my $referenceGenomeInfo        = "${pipelineDirName}/experiments/${experimentName}/pre_processing/inputInfo";

exit if (-d $referenceGenomeInfo);


#CREATE THE REFERENCE GENOME DIR
createReferenceGenomeDir();

#JUST PREPARING
if (($gtfName ne 'none') && ($fastaName ne 'none')){
 # (system "ln -s $gtfName ${referenceGenomeInfo}/___ex___.gtf")     == 0 or die "Error[createInput.pl]! Cannot execute:\nln -s $gtfName ${referenceGenomeInfo}/___ex___.gtf\n$!\n";
 # (system "ln -s $fastaName ${referenceGenomeInfo}/___query___.fa") == 0 or die "Error[createInput.pl]! Cannot execute:\nln -s $fastaName ${referenceGenomeInfo}/___query___.fa\n$!\n";
    (system "cp $gtfName ${referenceGenomeInfo}/___ex___.gtf")     == 0 or die "Error[createInput.pl]! Cannot execute:\ncp $gtfName ${referenceGenomeInfo}/___ex___.gtf\n$!\n";
    (system "cp $fastaName ${referenceGenomeInfo}/___query___.fa") == 0 or die "Error[createInput.pl]! Cannot execute:\ncp $fastaName ${referenceGenomeInfo}/___query___.fa\n$!\n";
  %infoGtf = readingGTFexons($gtfName);
}


#GTF2FASTA
if (($gtfName ne 'none') && ($fastaName eq 'none')){
  print STDERR "#EXTRACTING THE FASTA SEQUENCES FROM GTF COORDINATES\n";
  open (FA , ">${referenceGenomeInfo}/___query___.fa") or die "Error[createInput.pl]! cannot create ${referenceGenomeInfo}/___query___.fa $!\n";
  %infoGtf = readingGTFexons($gtfName);   #   copyExonGtf(); exit;
  foreach my $transcript (keys %infoGtf){
    my $string;
    my $fileChr = "${referenceGenomeInfo}/allGenomeInfo/*/chr/" . $infoGtf{$transcript}[0]{'chr'};
    #sequence
    @{$infoGtf{$transcript}} = sort {$a->{"start"}  <=>  $b->{"start"}} @{$infoGtf{$transcript}};
    foreach my $exon (0..$#{$infoGtf{$transcript}}){
      my $ex_start = $infoGtf{$transcript}[$exon]{'start'};
      my $ex_end   = $infoGtf{$transcript}[$exon]{'end'};
      $cmd = "$chr_subseqName $fileChr $ex_start $ex_end";
      $string .= `$cmd`; die ("Error[createInput.pl]! chr_subseq didn't work for the command $cmd $!\n") if ($?);
      chomp $string;
    }

    #strand
    if ($infoGtf{$transcript}[0]{'strand'} eq '-'){
      $string = getReverseComplement($string);
    }
    elsif ($infoGtf{$transcript}[0]{'strand'} eq '+'){}
    else {
      die "Error[createInput.pl]! the strand information in the gtf file must be either + or -. if you have \".\" you should either remove these lines, either replace it with + or -, or set \"-unstrandedGTF on\"  . Please fix it\n";
    }
    #else {unstrandedGTFfunction($transcript); next;}

    $string = string2FASTA($string);
    print FA ">$transcript\n$string\n";
  }
  close FA;
  copyExonGtf();
}

#FASTA2GTF
if (($gtfName eq 'none') && ($fastaName ne 'none')){
   print STDERR "#GENERATING THE QUERY GTF ANNOTATIONs ON THE REFERENCE GENOME\n";
   #(system "ln -s $fastaName ${referenceGenomeInfo}/___query___.fa") == 0 or die "Error[createInput.pl]! Cannot link the $fastaName file\n$!\n" ;
   (system "cp $fastaName ${referenceGenomeInfo}/___query___.fa") == 0 or die "Error[createInput.pl]! Cannot copy the $fastaName file\n$!\n" ;

   #blasting
   $cmd = "${pipelineDirName}/scripts/blastSearch.pl -query $fastaName -blast $blastName -blast_strategy $strategyName -splitGenome $splitGenome  -cluster $clusterName -experiment inputInfo  -pipeline_dir $referenceGenomeInfo";
   (system "$cmd") == 0 or die "Error[createInput.pl]! Cannot execute:\n$cmd\n$!\n";
   #exonerating
   cleaning_done();
   $cmd = "${pipelineDirName}/scripts/exonerateRemapping.pl  -mf2 ${referenceGenomeInfo}/experiments/inputInfo/BLAST_OUT/* -query $fastaName -targetGenomeFolder ${referenceGenomeInfo}/allGenomeInfo/*/chr  -exonerate_lines_mode $exonerate_lines_mode -exonerate_success_mode  1 -cluster $clusterName    -experiment inputInfo  -pipeline_dir $referenceGenomeInfo";
   (system "$cmd") == 0 or die "Error[createInput.pl]! Cannot execute:\n$cmd\n$!\n";
   waiting() if ($clusterName eq 'on');
   opendir (D,"${referenceGenomeInfo}/experiments/inputInfo/EXONERATE_OUT/");
   my @exonerateOutDir = readdir (D);
   closedir D;
   my $newGtf;
   foreach my $exonerateOutFile (@exonerateOutDir){
     if ($exonerateOutFile=~/\.gtf$/){
       $newGtf = "${referenceGenomeInfo}/experiments/inputInfo/EXONERATE_OUT/" . $exonerateOutFile;
     }
   }
   (system "mv $newGtf ${referenceGenomeInfo}/___ex___.gtf") == 0 or die "Error[createInput.pl]! Cannot move the $newGtf file\n$!\n";
   %infoGtf = readingGTFexons("${referenceGenomeInfo}/___ex___.gtf");
   cleaning_done();
 }

#EXONGTF2TRANSCRIPTGTF
exon2transcriptGtf();
#TRANSCRIPTGTF2GENOMICSIZE
transcriptGtf2genomicSize();


#UPDATE THE RNAmapping_pipeline.txt FILE
update_RNAmapping_pipeline();

















###
#FUNCTIONS
sub update_RNAmapping_pipeline {
  my $rnaMappingPipelineFile  = "${pipelineDirName}/experiments/${experimentName}/RNAmapping_pipeline.txt";
  my $tmpPipelineFile = fileNameGenerator("${rnaMappingPipelineFile}_tmp");
  open (NEW,">$tmpPipelineFile") or die "Error[createInput.pl]! Cannot create the $tmpPipelineFile file\n$!\n";
  open (RMP,"<$rnaMappingPipelineFile") or die "Error[createInput.pl]! Cannot read the $rnaMappingPipelineFile file\n$!\n";
  foreach my $line (<RMP>){
    chomp $line;
    if ($line=~/^QUERYFILE=(.+\s*)$/){
      $line = "QUERYFILE=${referenceGenomeInfo}/___query___.fa";
    }
    if ($line=~/^QUERYGTF=(.+\s*)$/){
      $line = "QUERYGTF=${referenceGenomeInfo}/___ex___.gtf";
    }
    if ($line=~/^PREPROCESSING=(.+\s*)$/){
      $line = "PREPROCESSING=off";
    }
    print NEW "$line\n";
  }
  close RMP;
  close NEW;
  (system "mv $tmpPipelineFile $rnaMappingPipelineFile") == 0 or die "Error[createInput.pl]! Cannot move $tmpPipelineFile to $rnaMappingPipelineFile\n$!\n";
}

sub cleaning_done {
  my $check_exonerate_done_ =  `ls $referenceGenomeInfo/experiments/inputInfo/CLUSTER_FILES/*___exonerate_done___ 2> /dev/null | wc -l`;
  chomp $check_exonerate_done_;
  if ($check_exonerate_done_ > 0){
    system "rm $referenceGenomeInfo/experiments/inputInfo/CLUSTER_FILES/*___exonerate_done___";
  }
}

sub help_message {
my $helpMessage = "\nNAME
createInput.pl - Input the queries in FASTA format (-fasta), exon gtf format (-gtf), or both. Input the reference genome as multi-FASTA file (-reference_genome). The output is stored in \"pre_processing/inputInfo\" experiment directory. The script will create the FASTA if not provided, and will create the gtf if not provided. These files will pass through the filteing.pl filters

SYNOPSIS
createInput.pl -experiment -pipeline_dir [-gtf -fasta -reference_genome -exonerate_mode -xdformat -blast_strategy -blast -cluster -chr_subseq ]

DESCRIPTION
   * The goal is having both FASTA and .gtf file to run the filtering and then the RNAmapping pipeline
   * If the user provides the exon .gtf file the FASTA sequences will be extracted from the reference genome unsing \"chr_subseq\"
   * If the user provides the FASTA file the script will call the RNAmapping pipeline and try to guess the .gtf file relying on the reference genome.
     This ammount in creating the the \"pre_processing/inputInfo\" RNAmapping pipeline working framework, with:
               -the reference_genome being formatted and stored in the subfolder \"allGenomeInfo\"
               -the blast and exonerate output stored in the experiment subfolder
   * The scripts will return also the transcript .gtf file and will replace the \"CONFIG/exonerateExtensionFile\" file with the genomic size of the queries

OPTIONS
   * -fasta                  <Default: none>
   * -gtf                    <Default: none>
   * -reference_genome       <Default: none>
   * -xdformat               <Default: xdformat>
   * -blast_strategy         <Default: wublastn_opt>
   * -blast                  <Default: wublastn>
   * -cluster                <Default: off>
   * -chr_subseq             <Default: chr_subseq>
   * -exonerate_success_mode <Default: 1>
   * -exonerate_lines_mode   <Default: 1000>

TROUBLESHOOTING
* The script will not be executed in the case \"pre_processing/inputInfo\" directory alredy exist. If you need to rerun createInput.pl on the same experiment you must first erease the \"pre_processing/inputInfo\" directory
";
print "$helpMessage\n\n\n";
exit;
}

sub options {
  my ($splitGenome , $unstrandedGTFname , $link2referenceGenomeName , $xdformatName , $gtfName , $fastaName , $referenceGenomeName , $chr_subseqName , $clusterName  , $strategyName , $blastName , $exonerate_success_mode , $exonerate_lines_mode , $experimentName , $pipelineDirName );
  my $spyUnstrandedGTF          = 1;
  my $spylink2referenceGenome   = 1;
  my $spyXdformat               = 1;
  my $spyPipelineDir            = 1;
  my $spyExperiment             = 1;
  my $spyGtfName                = 1;
  my $spyFastaName              = 1;
  my $spyReferenceGenomeName    = 1;
  my $spyChr_subseq             = 1;
  my $spyCluster                = 1;
  my $spyStrategy               = 1;
  my $spyBlast                  = 1;
  my $spyExonerate_lines_mode   = 1;
  my $spyExonerate_success_mode = 1;
  my $spySplitGenome            = 1;


  foreach my $field (0..$#ARGV){
    if ($ARGV[$field] eq '-splitGenome'){
      $splitGenome = $ARGV[1+$field];
      $spySplitGenome = 2;
      next;
    }
    if ($spySplitGenome == 2){
      $spySplitGenome = 3;
      next;
    }
    if ($ARGV[$field] eq '-unstrandedGTF'){
      $unstrandedGTFname = $ARGV[1+$field];
      $spyUnstrandedGTF = 2;
      next;
    }
    if ($spyUnstrandedGTF == 2){
      $spyUnstrandedGTF = 3;
      next;
    }
    if ($ARGV[$field] eq '-link2referenceGenome'){
      $link2referenceGenomeName = $ARGV[1+$field];
      $spylink2referenceGenome = 2;
      next;
    }
    if ($spylink2referenceGenome == 2){
      $spylink2referenceGenome = 3;
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
    if ($ARGV[$field] eq '-experiment'){
	$experimentName = $ARGV[1+$field];
	$spyExperiment = 2;
	next;
    }
    if ($spyExperiment == 2){
	$spyExperiment = 3;
	next;
    }
    if ($ARGV[$field] eq '-exonerate_success_mode'){
	$exonerate_success_mode = $ARGV[1+$field];
	$spyExonerate_success_mode = 2;
	next;
    }
    if ($spyExonerate_success_mode == 2){
      $spyExonerate_success_mode = 3;
      next;
    }
    if ($ARGV[$field] eq '-exonerate_lines_mode'){
	$exonerate_lines_mode = $ARGV[1+$field];
	$spyExonerate_lines_mode = 2;
	next;
    }
    if ($spyExonerate_lines_mode == 2){
      $spyExonerate_lines_mode = 3;
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
    if ($ARGV[$field] eq '-cluster'){
	$clusterName = $ARGV[1+$field];
	$spyCluster = 2;
	next;
    }
    if ($spyCluster == 2){
	$spyCluster = 3;
	next;
    }
    if ($ARGV[$field] eq '-chr_subseq'){
      $chr_subseqName = $ARGV[1+$field];
      $spyChr_subseq = 2;
      next;
    }
    if ($spyChr_subseq == 2){
      $spyChr_subseq = 3;
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
    if ($ARGV[$field] eq '-gtf'){
      $gtfName = $ARGV[1+$field];
      $spyGtfName = 2;
      next;
    }
    if ($spyGtfName == 2){
      $spyGtfName = 3;
      next;
    }
    if ($ARGV[$field] eq '-fasta'){
      $fastaName = $ARGV[1+$field];
      $spyFastaName = 2;
      next;
    }
    if ($spyFastaName == 2){
      $spyFastaName  = 3;
      next;
    }
    if ($ARGV[$field] eq '-reference_genome'){
	$referenceGenomeName = $ARGV[1+$field];
	$spyReferenceGenomeName = 2;
	next;
    }
    if ($spyReferenceGenomeName == 2){
	$spyReferenceGenomeName = 3;
	next;
    }
  }

  die "Error[createInput.pl]! You must provide the -pipeline_dir parameter.\n"                                    if ($spyPipelineDir  != 3);
  die "Error[createInput.pl]! You must provide the -experiment parameter.\n"                                      if ($spyExperiment  != 3);
  die "Error[createInput.pl]! The -cluster parameter can be either \"on\" or \"off\"\n"                           if ((defined $clusterName) && ($clusterName ne 'on')  && ($clusterName ne 'off'));
  die "Error[blastSearch.pl]! -splitGenome parameter can be either \'yes\' or \'no\'\n"                           if ((defined $splitGenome) && ($splitGenome ne 'yes') && ($splitGenome ne 'no'));
  die "Error[blastSearch.pl]! The -splitGenome function is meant to run on the cluster, please set -cluster on\n" if ((defined $splitGenome) && ($splitGenome eq 'yes') && ($clusterName ne 'on'));

  $unstrandedGTFname        = 'off'          if (! defined $unstrandedGTFname);
  $link2referenceGenomeName = 'off'          if (! defined $link2referenceGenomeName);
  $fastaName                = 'none'         if (! defined $fastaName);
  $gtfName                  = 'none'         if (! defined $gtfName);
  $chr_subseqName           = 'chr_subseq'   if (! defined $chr_subseqName);
  $clusterName              = 'off'          if (! defined $clusterName);
  $strategyName             = 'wublastn_opt' if (! defined $strategyName);
  $blastName                = 'wu-blastn'    if (! defined $blastName);
  $xdformatName             = 'xdformat'     if (! defined $xdformatName);
  $splitGenome              = 'no'           if (! defined $splitGenome);

  die "Error[createInput.pl]! You should specify at least a FASTA or .gtf file\n" if (($fastaName eq 'none') && ($gtfName eq 'none'));
  die "Error[createInput.pl]! You must provide the -reference_genome parameter.\n" if (($spyReferenceGenomeName != 3) && (($fastaName eq 'none') || ($gtfName eq 'none')));
  if ((defined $exonerate_lines_mode) && ($exonerate_lines_mode ne "exhaustive") && ($exonerate_lines_mode =~/\D/)){
      die "Error[exonerateRemapping.pl]! The supported exonerate_lines_mode modes are either \"exhaustive\" either an integer number of iteration\n";
  }
  if ((defined $exonerate_success_mode) && ($exonerate_success_mode ne "exhaustive") && ($exonerate_success_mode ne "ortholog") && ($exonerate_success_mode =~/\D/)){
      die "Error[exonerateRemapping.pl]! The supported exonerate_success_mode modes are either \"exhaustive\" either \"ortholog\" either  an integer number of iteration\n";
  }
  $exonerate_lines_mode   = 1000  if (! defined $exonerate_lines_mode);
  $exonerate_success_mode = 1     if (! defined $exonerate_success_mode);

  return ($splitGenome , $unstrandedGTFname , $link2referenceGenomeName , $xdformatName , $gtfName , $fastaName , $referenceGenomeName , $chr_subseqName , $clusterName  , $strategyName , $blastName , $exonerate_success_mode , $exonerate_lines_mode  , $experimentName , $pipelineDirName );
}

sub acceptedVariableSpace {
  my %space = ('-splitGenome' => 1 , '-unstrandedGTF' => 1 , '-link2referenceGenome' => 1 , '-xdformat' => 1 , '-experiment' => 1 ,  '-blast_strategy' => 1 ,  '-blast' => 1 ,  '-cluster' => 1 ,  '-chr_subseq' => 1 ,  '-pipeline_dir' => 1 ,  '-gtf' => 1 ,  '-fasta' => 1 ,  '-reference_genome' => 1 , '-exonerate_lines_mode' => 1  , '-exonerate_success_mode' => 1  );
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[createInput.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}



#read the hash returned by the function readingGTFexons() and returns 3 hash references, one having the gene sizes, one having transcript sizes and one having exon sizes
sub gtfExons2sizes {
  my %infoExons          =  @_;
  my (%infoExonsPerGene , %transcriptSizes) ;
  my (@transcriptNames, @geneNames);

  foreach my $transcript (keys %infoExons){
    foreach my $exon (0..$#{$infoExons{$transcript}}){
      my $chr      = $infoExons{$transcript}[$exon]->{'chr'};
      my $start    = $infoExons{$transcript}[$exon]->{'start'};
      my $end      = $infoExons{$transcript}[$exon]->{'end'};
      my $gene     = $infoExons{$transcript}[$exon]->{'gene_id'};
      my $rna      = $infoExons{$transcript}[$exon]->{'transcript_id'};
      my $exonSize = abs ($end - $start);
      $infoExons{$transcript}[$exon]->{'size'} = $exonSize;

      push (@transcriptNames , $transcript) if (! defined $transcriptSizes{$transcript});
      push (@geneNames, $gene)       if (! defined $infoExonsPerGene{$gene});
      $transcriptSizes{$transcript} += $exonSize;

      my %hash = ("chr" =>  $chr, "start" => $start, "end" => $end);
      push (@{$infoExonsPerGene{$gene}}, \%hash);
    }
  }

  #TAKING GENE SIZES
  my %geneSizes;
  foreach my $gene (keys %infoExonsPerGene){
    my $minExonStart = 1000000000000000000000000000000000000000000000000000000000000000000000000000;
    my $maxExonEnd   = 0;
    my $chrKeeper;
    foreach my $index (0..$#{$infoExonsPerGene{$gene}}){
      $minExonStart = min ($minExonStart, $infoExonsPerGene{$gene}->[$index]->{"start"});
      $maxExonEnd   = max ($maxExonEnd, $infoExonsPerGene{$gene}->[$index]->{"end"});
      $chrKeeper    = $infoExonsPerGene{$gene}->[$index]->{"chr"};
    }
    my $geneSize = abs ($maxExonEnd - $minExonStart);
    $geneSizes{$gene} = $geneSize;
  }

  return (\%geneSizes , \%transcriptSizes , \%infoExons);

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


sub getReverseComplement {
  my ($sequence) = @_;
  $sequence = reverse $sequence;
  my $complementedSequence = "";
  for (my $key = 0; $key < length($sequence); $key++) {
    my $char = substr $sequence, $key, 1;
    if    ($char eq "A") {$complementedSequence .= "T";}
    elsif ($char eq "C") {$complementedSequence .= "G";}
    elsif ($char eq "G") {$complementedSequence .= "C";}
    elsif ($char eq "T") {$complementedSequence .= "A";}
    elsif ($char eq "U") {$complementedSequence .= "A";}
    elsif ($char eq "a") {$complementedSequence .= "t";}
    elsif ($char eq "c") {$complementedSequence .= "g";}
    elsif ($char eq "g") {$complementedSequence .= "c";}
    elsif ($char eq "t") {$complementedSequence .= "a";}
    elsif ($char eq "u") {$complementedSequence .= "a";}
    elsif ($char eq "(") {$complementedSequence .= ")";}
    elsif ($char eq ")") {$complementedSequence .= "(";}
    elsif ($char eq "<") {$complementedSequence .= ">";}
    elsif ($char eq ">") {$complementedSequence .= "<";}
    else {
        $complementedSequence .= $char;
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

sub waiting {
    my $total = 1;
    my $time  = 10;
    my $done  = 0;
    while ($done < $total ){
	print STDERR "#Waiting the job to exit the cluster... " . $done . " already finished\n";
	system "sleep $time";
	$time += 20;
	$done = `ls ${referenceGenomeInfo}/experiments/inputInfo/CLUSTER_FILES/*___exonerate_done___ 2> /dev/null | wc -l`;
	chomp $done;
    }
    print STDERR "#All jobs exited the cluster\n";
}

sub gtfExons2Transcripts{
  my (%infoGtf) = @_;
  my %infoGtfTranscript;
  foreach my $rna (keys %infoGtf){
    my $minExonStart = 1000000000000000000000000000000000000000000000000000000000000000000000000000;
    my $maxExonEnd   = 0;
    my ($chr , $group , $frame , $strand , $score , $feature , $source , $gene_id);
    foreach my $index (0..$#{$infoGtf{$rna}}){
      $minExonStart = min2 ($minExonStart, $infoGtf{$rna}->[$index]->{"start"});
      $maxExonEnd   = max2 ($maxExonEnd, $infoGtf{$rna}->[$index]->{"end"});

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

sub exon2transcriptGtf {
  open (O,">${referenceGenomeInfo}/___tx___.gtf") or die "Error[createInput]! cannot create the ${referenceGenomeInfo}/___tx___.gtf transcript gtf output\n";
  my %infoTranscripts = gtfExons2Transcripts(%infoGtf);
  foreach my $id (keys %infoGtf){
    my $outLine;
    $outLine .= $infoTranscripts{$id}->{'chr'}     . "\t";
    $outLine .= $infoTranscripts{$id}->{'source'}  . "\t";
    $outLine .= "transcript" . "\t";
    $outLine .= $infoTranscripts{$id}->{'start'}   . "\t";
    $outLine .= $infoTranscripts{$id}->{'end'}     . "\t";
    $outLine .= $infoTranscripts{$id}->{'score'}   . "\t";
    $outLine .= $infoTranscripts{$id}->{'strand'}  . "\t";
    $outLine .= $infoTranscripts{$id}->{'frame'}   . "\t";
    $outLine .= $infoTranscripts{$id}->{'group'}   . "\t";

    print O "$outLine\n";
  }
  close O;
}



sub transcriptGtf2genomicSize {
  open (O,"<${referenceGenomeInfo}/___tx___.gtf") or die "Error[createInput]! cannot read the ${referenceGenomeInfo}/___tx___.gtf transcript gtf output\n";
  open (C, ">${pipelineDirName}/experiments/${experimentName}/CONFIG/exonerateExtensionFile") or die "Error[createInput]! cannot create the transcript size output\n";
  foreach my $line (<O>){
    my ($tx_id , $size);
    if ($line =~/transcript_id \"([^\"]+)\"/){
      $tx_id = $1;
    }
    else{
      die "Error[createInput]! ${pipelineDirName}/QUERY/___tx___.gtf does not have a transcript_id field\n";
    }
    if ($line =~/^\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+/){
      $size = abs ($2 - $1);
    }
    print C "$tx_id $size\n";
  }
  close O;
  close C;
}


sub createReferenceGenomeDir {
  system "mkdir $referenceGenomeInfo"                                                                                unless (-d $referenceGenomeInfo);
  system "mkdir $referenceGenomeInfo/experiments"                                                                    unless (-d "${referenceGenomeInfo}/experiments");
  system "mkdir $referenceGenomeInfo/allGenomeInfo"                                                                  unless (-d "${referenceGenomeInfo}/allGenomeInfo");
  system "mkdir $referenceGenomeInfo/experiments/inputInfo"                                                          unless (-d "${referenceGenomeInfo}/experiments/inputInfo");
  system "cp -r ${pipelineDirName}/experiments/${experimentName}/CONFIG  $referenceGenomeInfo/experiments/inputInfo" unless (-d "${referenceGenomeInfo}/experiments/inputInfo/CONFIG");
  system "mkdir $referenceGenomeInfo/experiments/inputInfo/CLUSTER_FILES"                                            unless (-d "${referenceGenomeInfo}/experiments/inputInfo/CLUSTER_FILES");
  system "mkdir $referenceGenomeInfo/experiments/inputInfo/STDERR"                                                   unless (-d "${referenceGenomeInfo}/experiments/inputInfo/STDERR");
  system "mkdir $referenceGenomeInfo/experiments/inputInfo/BLAST_OUT"                                                unless (-d "${referenceGenomeInfo}/experiments/inputInfo/BLAST_OUT");
  system "mkdir $referenceGenomeInfo/experiments/inputInfo/EXONERATE_OUT"                                            unless (-d "${referenceGenomeInfo}/experiments/inputInfo/EXONERATE_OUT");

  (system "cp ${pipelineDirName}/experiments/${experimentName}/RNAmapping_pipeline.txt $referenceGenomeInfo/experiments/inputInfo/")==0 or die "Error[createInput]! cannot copy the RNAmapping_pipeline.txt in the pre_processing/experiments/inputInfo folder\n$!\n";

  if ($link2referenceGenomeName eq 'off'){
    #my $referenceGenomeBaseName = `basename $referenceGenomeName`; if ($?){die "Error[createInput.pl]! Cannot run basename $referenceGenomeName \n$!\n";}
    my $genomePointer = fileNameGenerator ("${referenceGenomeInfo}/___referenceGenomePointer___");
    open (GP,">$genomePointer") or die "Error[createInput.pl]! cannot create $genomePointer $!\n";
    #print GP "$referenceGenomeName $referenceGenomeBaseName";
    print GP "$referenceGenomeName referenceGenome";
    close GP;
    $cmd = "${pipelineDirName}/scripts/createDatabase.pl  -genomes $genomePointer -xdformat $xdformatName -splitGenome $splitGenome  -strategy $strategyName -pipeline_dir $referenceGenomeInfo";
    (system "$cmd") == 0 or die "Error[createInput.pl]! Cannot run $cmd\n$!\n";
    system "rm $genomePointer";
  }
  else{
     if (-e "${pipelineDirName}/allGenomeInfo/$link2referenceGenomeName"){
       if (-l "${pipelineDirName}/allGenomeInfo/$link2referenceGenomeName"){
	 my $linkLine = `ls -l ${pipelineDirName}/allGenomeInfo/$link2referenceGenomeName`; if ($?){die "Error[createInput.pl]! failed with:\nls -l ${pipelineDirName}/allGenomeInfo/$link2referenceGenomeName\n$!\n";}
	 chomp $linkLine;
	 if ($linkLine=~/-> (\S+)$/){
	   my $linkAddress = $1;
	   (system "ln -s $linkAddress ${referenceGenomeInfo}/allGenomeInfo/${link2referenceGenomeName}") == 0 or die "Error[createInput.pl]! cannot create the symbolic link. The command that failed was:\nln -s $linkAddress ${referenceGenomeInfo}/allGenomeInfo/${link2referenceGenomeName}\n$!\n";
	 }
       }
       else{
     	(system "ln -s ${pipelineDirName}/allGenomeInfo/$link2referenceGenomeName ${referenceGenomeInfo}/allGenomeInfo/${link2referenceGenomeName}") == 0 or die "Error[createInput.pl]! cannot create the symbolic link. The command that failed was:\nln -s ${pipelineDirName}/allGenomeInfo/${link2referenceGenomeName} ${referenceGenomeInfo}/allGenomeInfo/${link2referenceGenomeName}\n$!\n";
       return;
       }
     }
     else{
       die "Error[createInput.pl]! cannot create the symbolic link. $link2referenceGenomeName does not exist\n";
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


sub readingGTFexons {
  my ($GTFfile) = @_;
  open (IN,"<$GTFfile") or die "Error[createInput.pl]! Cannot open the GTFfile:\n$GTFfile \n$!";
  my %infoGtf;
  foreach my $line (<IN>){
	chomp $line;
	my ($chr , $group , $frame , $strand , $score , $end , $start , $feature , $source , $gene_id , $transcript_id);
	if ($line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.*)$/){
	    $chr     = $1;
	    $source  = $2;
	    $feature = $3;
	    $start   = min2 ($4, $5);
	    $end     = max2 ($4, $5);
	    $score   = $6;
	    $strand  = $7;
	    $frame   = $8;
	    $group   = $9;
	    if ($group =~/gene_id \"([^\"]+)\"/){
		$gene_id = $1;
	    }
	    else {die "Error[createInput.pl]! gtf file in wrong format. Impossible to find the gene_id in line:\n$line\n when it is mandatory for this format\n";}
	    if ($group =~/transcript_id \"([^\"]+)\"/){
		$transcript_id = $1;
	    }
	    else {die "Error[createInput.pl]! gtf file in wrong format. Impossible to find the transcript_id in line:\n$line\n when it is mandatory for this format\n";}

	    next unless ($feature eq 'exon');


	    if (($unstrandedGTFname eq 'on') && ($strand eq '.')){
	      $group =~s/transcript_id \"[^\"]+\"/transcript_id \"${transcript_id}_noStrandPlus\"/;
	      $group =~s/gene_id \"[^\"]+\"/gene_id \"${gene_id}_noStrandPlus\"/;
	      my  %hash1 = (
			   "chr"           => $chr,
			   "source"        => $source,
			   "feature"       => $feature,
			   "start"         => $start,
			   "end"           => $end,
			   "score"         => $score,
			   "strand"        => '+',
			   "frame"         => $frame,
			   "group"         => $group,
			   "gene_id"       => $gene_id,
			  );
	      push (@{$infoGtf{"${transcript_id}_noStrandPlus"}}, \%hash1);
	      $group =~s/transcript_id \"[^\"]+\"/transcript_id \"${transcript_id}_noStrandMinus\"/;
	      $group =~s/gene_id \"[^\"]+\"/gene_id \"${gene_id}_noStrandMinus\"/;
	      my  %hash2 = (
			   "chr"           => $chr,
			   "source"        => $source,
			   "feature"       => $feature,
			   "start"         => $start,
			   "end"           => $end,
			   "score"         => $score,
			   "strand"        => '-',
			   "frame"         => $frame,
			   "group"         => $group,
			   "gene_id"       => $gene_id,
			  );
	      push (@{$infoGtf{"${transcript_id}_noStrandMinus"}}, \%hash2);

	    }
	    else{
	      my  %hash = (
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

	}
	else {die "Error[createInput.pl]! cannot read the line $line\n";}
      }
  close IN;
  return %infoGtf;
}
# sub unstrandedGTFfunction {
#   my ($transcript) = @_;
#   $infoGtf{"${transcript}_noStrandPlus"}  = $infoGtf{$transcript};
#   foreach my $exon (0..$#{$infoGtf{"$transcript"}}){
#     foreach my $field (keys %{$infoGtf{"$transcript"}[$exon]}){
#       my $value = $infoGtf{"$transcript"}[$exon]{$field};
#       $infoGtf{"${transcript}_noStrandMinus"}[$exon]{$field} = $value;
#     }
#   }
#   delete $infoGtf{$transcript};

#   foreach my $exon (0..$#{$infoGtf{"${transcript}_noStrandPlus"}}){
#     $infoGtf{"${transcript}_noStrandPlus"}[$exon]->{'strand'} = '+';
#   }
#   foreach my $exon (0..$#{$infoGtf{"${transcript}_noStrandMinus"}}){
#     $infoGtf{"${transcript}_noStrandMinus"}[$exon]->{'strand'} = '-';
#   }
# }

sub copyExonGtf {   #foreach my $id (keys %infoGtf){print "\'$id\'\n";}
  open (EX , ">${referenceGenomeInfo}/___ex___.gtf") or die "Error[createInput.pl]! cannot create ${referenceGenomeInfo}/___ex___.gtf\n$!\n";
  foreach my $id (keys %infoGtf){
    foreach my $exon (0..$#{$infoGtf{$id}}){
      my $outLine;
      $outLine .= $infoGtf{$id}[$exon]->{'chr'}     . "\t";
      $outLine .= $infoGtf{$id}[$exon]->{'source'}  . "\t";
      $outLine .= $infoGtf{$id}[$exon]->{'feature'} . "\t";
      $outLine .= $infoGtf{$id}[$exon]->{'start'}   . "\t";
      $outLine .= $infoGtf{$id}[$exon]->{'end'}     . "\t";
      $outLine .= $infoGtf{$id}[$exon]->{'score'}   . "\t";
      $outLine .= $infoGtf{$id}[$exon]->{'strand'}  . "\t";
      $outLine .= $infoGtf{$id}[$exon]->{'frame'}   . "\t";
      $outLine .= $infoGtf{$id}[$exon]->{'group'}   . "\t";
      print EX "$outLine\n";
    }
  }
  close EX;
}
