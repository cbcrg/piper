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

#  Usage: motevo input_sequences param_file wm_file
#TAKE OPTIONS
acceptedVariableSpace();
my ($type , $repeatThreshold , $motevo , $motevo_tree , $motevo_wm_file , $motevo_runUFE , $motevo_reference , $t_coffee  ,  $experimentName , $pipelineDirName) = options();
my $outDir        = "${pipelineDirName}/experiments/${experimentName}/qualityCheck/alignPromoters/";
my $promotersDir  = "${pipelineDirName}/experiments/${experimentName}/qualityCheck/extractedFeatues/promoters_${type}_rep${repeatThreshold}/";
die "Error[promoters_check.pl]! To run this script it is needed the folder\n$promotersDir\nYou can have it by running extractFeatures\n" unless (-d $promotersDir);



#exon gtf 2 promoter gtf
(system "mkdir -p  $outDir") == 0 or die "Error[promoters_check.pl]! cannot execute\nmkdir -p  $outDir\n$!\n";
printLog($outDir);
opendir (E,"$promotersDir") or die "Error[promoters_check.pl]! Cannot open $promotersDir\n$!\n";
while (my $prFastaFile = readdir(E)){
  next if (($prFastaFile eq '.') || ($prFastaFile eq '..') || ($prFastaFile!~/\.fa$/));
  my $prAlnFile = $prFastaFile;
  $prAlnFile =~s/\.fa$/\.aln/;
  my $cmd = "$t_coffee -mode procoffee -in ${promotersDir}/${prFastaFile} -quiet -outfile ${outDir}/${prAlnFile}";
  (system "$cmd") == 0 or die "Error[promoters_check.pl]! Cannot run\n$cmd\n$!\n";
}
closedir E;
(system "rm ${outDir}/*.html") == 0 or die "Error[promoters_check.pl]! Cannot execute\nrm ${outDir}/*.html\n$!\n";
#convert to fasta_aln format
opendir (O,"$outDir") or die "Error[promoters_check.pl]! Cannot open $outDir\n$!\n";
while (my $f = readdir(O)){
  chomp $f;
  next if (($f eq '.') || ($f eq '..') || ($f eq 'log.txt'));
  my $cmd = "$t_coffee -other_pg seq_reformat -in ${outDir}/$f -output fasta_aln > ${outDir}/${f}_2";
  (system "$cmd") == 0 or die "Error[promoters_check.pl]! Cannot execute\n$cmd\n$!\n";
  (system "mv ${outDir}/${f}_2 ${outDir}/${f}") == 0 or die "Error[promoters_check.pl]! Cannot execute\nmv ${outDir}/${f}_2 ${outDir}/${f}\n$!\n";
}
closedir O;




#TFBS analysis
my $TFBS_analysis = TFBS_analysis_check();
my (%bg , $tree , $parameterFile);
my $outMotevo        = "${pipelineDirName}/experiments/${experimentName}/qualityCheck/MotEvo/";
(system "mkdir -p  $outMotevo") == 0 or die "Error[promoters_check.pl]! cannot execute\nmkdir -p  $outMotevo\n$!\n";
if ($TFBS_analysis == 1){

  #read tree
  $tree= readTree();

  #take background frequences
  %bg = takeBGfreq();

  #UFE
  my $UFEcmd = "$motevo_runUFE $motevo_tree " . $bg{'A'} . " " . $bg{'C'} . " ". $bg{'G'} . " " .   $bg{'T'} . " > ${outMotevo}/UFE";
  (system "$UFEcmd") == 0 or die "Error[promoters_check.pl]! Cannot execute\n$UFEcmd\n$!\n" ;                                                        #trick
  my @species = takeSpecies("${outMotevo}/UFE");

  #make parameter File
  $parameterFile = makeParameterFile();

  opendir (O,"$outDir") or die "Error[promoters_check.pl]! Cannot open $outDir\n$!\n";
  while (my $f = readdir(O)){
      next if (($f eq '.') || ($f eq '..'));
      my $check = trimUnusedSpecies(\@species , $f);
      next if($check eq "___skipp___");
      my $cmd = "$motevo ${outMotevo}/$f  $parameterFile  $motevo_wm_file > ${outMotevo}/${f}.mot";
      (system "$cmd") == 0 or die "Error[promoters_check.pl]! Cannot execute\n$cmd\n$!\n";
      (system "mv priors_ma ${outMotevo}/${f}.priors.mot") == 0 or die "Error[promoters_check.pl]! Cannot execute:\nmv priors_ma ${outMotevo}/${f}.priors.mot\n$!\n"   ;
      (system "mv sites_ma ${outMotevo}/${f}.sites.mot") == 0 or die "Error[promoters_check.pl]! Cannot execute:\nmv sites_ma ${outMotevo}/${f}.sites.mot\n$!\n"  ;
   }
   closedir O;
}














#FUNCTIONS
sub options {
  my ($type , $repeatThreshold , $motevo , $motevo_tree_file , $motevo_wm_file , $motevo_runUFE , $motevo_reference , $t_coffee , $experimentName , $pipelineDirName);
  my $spyExperiment           = 1;
  my $spyPipelineDir          = 1;
  my $spyT_coffee             = 1;
  my $spyMotevo               = 1;
  my $spyMotevoTreeFile       = 1;
  my $spyMotevoWmFile         = 1;
  my $spyMotevoRunUFE         = 1;
  my $spyMotevoReference      = 1;
  my $spyRepeatThreshold      = 1;
  my $spyType                 = 1;

  foreach my $field (0..$#ARGV){
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
    if ($ARGV[$field] eq '-t_coffee'){
	$t_coffee = $ARGV[1+$field];
	$spyT_coffee = 2;
	next;
    }
    if ($spyT_coffee == 2){
	$spyT_coffee = 3;
	next;
    }
    if ($ARGV[$field] eq '-motevo'){
	$motevo = $ARGV[1+$field];
	$spyMotevo = 2;
	next;
    }
    if ($spyMotevo == 2){
	$spyMotevo = 3;
	next;
    }
    if ($ARGV[$field] eq '-motevo_reference'){
	$motevo_reference = $ARGV[1+$field];
	$spyMotevoReference = 2;
	next;
    }
    if ($spyMotevoReference == 2){
	$spyMotevoReference = 3;
	next;
    }
    if ($ARGV[$field] eq '-motevo_runUFE'){
	$motevo_runUFE = $ARGV[1+$field];
	$spyMotevoRunUFE = 2;
	next;
    }
    if ($spyMotevoRunUFE == 2){
	$spyMotevoRunUFE = 3;
	next;
    }
    if ($ARGV[$field] eq '-motevo_tree_file'){
	$motevo_tree_file = $ARGV[1+$field];
	$spyMotevoTreeFile = 2;
	next;
    }
    if ($spyMotevoTreeFile == 2){
	$spyMotevoTreeFile = 3;
	next;
    }
    if ($ARGV[$field] eq '-motevo_wm_file'){
        $motevo_wm_file = $ARGV[1+$field];
	$spyMotevoWmFile = 2;
	next;
    }
    if ($spyMotevoWmFile == 2){
	$spyMotevoWmFile = 3;
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
  die "Error[promoters_check.pl]! Please provide the -experiment parameter\n"   if (! defined $experimentName);
  die "Error[promoters_check.pl]! Please provide the -pipeline_dir parameter\n" if (! defined $pipelineDirName);
  die "Error[promoters_check.pl]! Please provide the -repeatThreshold parameter\n" if (! defined $repeatThreshold);
  die "Error[promoters_check.pl]! Please provide the -type parameter\n"            if (! defined $type);

  #defaults
  $t_coffee = 't_coffee' if (! defined $t_coffee);
  return ($type , $repeatThreshold , $motevo , $motevo_tree_file , $motevo_wm_file , $motevo_runUFE , $motevo_reference , $t_coffee  ,  $experimentName , $pipelineDirName);
}

sub acceptedVariableSpace {
  my %space = ('-type' => 1 , '-repeatThreshold' => 1 , '-motevo' => 1 , '-motevo_tree_file' =>1 , '-motevo_wm_file' => 1 , '-experiment' => 1 , '-pipeline_dir' => 1   , '-t_coffee' => 1 , '-motevo_runUFE' => 1 , '-motevo_reference' => 1);
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[promoters_check.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}

sub help_message {
my $helpMessage = "\nNAME
promoters_check.pl - it extracts the region upstream the mapped transcripts and check for promoter motifs

SYNOPSIS
promoters_check.pl -experiment -pipeline_dir [-promoters_size] [-motevo -motevo_tree_file -motevo_wm_file -motevo_runUFE -motevo_reference]

DESCRIPTION
   * WARNING1: if you run the Motevo prediction bear in mind that the UFE must be estimated on a tree where the maximum number of species recommended is 10 (otherwise it takes too much)
   * To get a tree of life containing just the species you want (max. 10) you can use the prune function of t_coffee, passing in -in the full tree of life (provided with pipeR)
   * and in -in2 a list of species name to pick up, in this format:
   *   >human
   *   >pig
   *   >cow
   * The command line is like
   *t_coffee -other_pg seq_reformat -action +tree_prune -in scripts/utility/pruned_tree_of_life -in2 speciesListFile
   * IMPORTANT!! Mind to remove the root of the tree if this is rooted. This is the last number to the right.
   * for instance, this pruned tree of life  ((human:0.00670,(mouse:0.08451,rat:0.09159):0.19777):0.02059,((cow:0.06180,pig:0.07900):0.02017,dog:0.05123):0.01167):0.02366;
   * should be manually changed in ((human:0.00670,(mouse:0.08451,rat:0.09159):0.19777):0.02059,((cow:0.06180,pig:0.07900):0.02017,dog:0.05123):0.01167);

   * The species of the MSA that will be used for motif prediction are the ones specified by the Newick tree.
   * MSA containing other species will be temporary re-written so that just the species specified in the motevo-tree are used

   * WARNING2: the specified species tree must have species name that can be found in the MSA headers
   * DEPENDENCY: To execute this check it is needed before to run successfully the pipeline up to the exonerate step

";
print "$helpMessage\n\n\n";
exit;
}

sub TFBS_analysis_check {
  if ((defined $motevo) && (defined $motevo_tree) && (defined $motevo_wm_file) && (defined $motevo_runUFE) && (defined $motevo_reference)){
    return 1;
  }
  if ((defined $motevo) || (defined $motevo_tree) || (defined $motevo_wm_file) || (defined $motevo_runUFE) || (defined $motevo_reference)){
    print "Error[promoters_check.pl]! To run motevo are needed at least 5 parameters. You can find the motevo and runUFE executables in motevo installation\nYou will find also the Weighted Matrices and the trees\nYou have to choose the reference species name\n";
    return 0;
  }
}
sub readTree {
  my $tree;
  open (T, "<$motevo_tree") or die "Error[promoters_check.pl] cannot open the tree file $motevo_tree\n$!\n";
  foreach my $l (<T>){
    chomp $l;
    next if ($l=~/^\s*$/);
    $tree = $l;
  }
  close T;
  return $tree;
}



sub takeSpecies {
  my ($UFE) = @_;
  my @species;
  open (UFE,"<$UFE") or die "Error[promoters_check.pl]! Cannot open $UFE\n$!\n";
  foreach my $l (<UFE>){
    chomp $l;
    if ($l=~/^>(\S+)/){
      my $s = $1;
      push (@species , $s);
    }
    else {
      last;
    }
  }
  close UFE;
  return @species;
}

sub trimUnusedSpecies {
  my ($species_p , $msa ) = @_;
  my @species = @{$species_p};
  my $spy   = 0;
  my $check = 0;
  my $msg = "___skipp___";
  open (F,"<${outDir}/$msa")    or die "Error[promoters_check.pl]! Cannot open ${outDir}/$msa\n$!\n";
  open (OU,">${outMotevo}/$msa") or die "Error[promoters_check.pl]! Cannot create ${outMotevo}/$msa\n$!\n";
  while (<F>){
    my $l = $_;
    chomp $l;
    if ($l=~/^>(\S+)/){
      my $id = $1;
      $spy   = 0;
      $check = 0;
      foreach my $s (@species){
	if ($id=~/$s/){
	  $spy = 1;
	  $check++;
	  if ($s eq $motevo_reference){
	    $spy   = 2;
	    $msg = "___ok___";
	  }
	}
      }
      die "Error[promoters_check.pl]! In ${outDir}/$msa there is this sequence\n$l\nthat match with more than one species in the tree\n$!\n" if ($check > 1);
      if ($spy == 2){
	print OU '>' . $l . "\n";
	next;
      }
    }
    if ($spy == 2){
      print OU $l . "\n";
    }
  }
  close F;

  if ($msg eq "___skipp___"){
    system "rm ${outMotevo}/$msa > /dev/null 2> /dev/null";
    return $msg;
  }


  open (F,"<${outDir}/$msa")    or die "Error[promoters_check.pl]! Cannot open ${outDir}/$msa\n$!\n";
  while (<F>){
    my $l = $_;
    chomp $l;
    if ($l=~/^>(\S+)/){
      my $id = $1;
      $spy = 0;
      foreach my $s (@species){
	if ($id=~/$s/){
	  $spy = 1;
	  if ($s eq $motevo_reference){
	    $spy   = 2;
	  }
	}
      }
      if ($spy == 1){
	  print OU $l . "\n";
	  next;
      }
    }
    if ($spy == 1){
      print OU $l . "\n";
      next;
    }
  }
  close F;
  close OU;
  return $msg;
}

sub makeParameterFile {
  open (P,">${outMotevo}/parameters") or die "Error[promoters_check.pl]! Cannot open ${outMotevo}/parameters\n$!\n";
  print P "refspecies $motevo_reference\n";
  print P "TREE $tree\n\n";
  print P "Mode TFBS\n";
  print P "EMprior 0\n\n";
  print P "UFEwmprior 200.0\n";
  print P "UFEwmfile ${outMotevo}/UFE\n";
  print P "UFEwmlen 8\n";
  print P "UFEprint 1\n\n";
  print P "markovorderBG 0\n";
  print P "bgprior 0.98\n";
  print P "bg A " . $bg{'A'} . "\n";
  print P "bg T " . $bg{'T'} . "\n";
  print P "bg G " . $bg{'G'} . "\n";
  print P "bg C " . $bg{'C'} . "\n";
  print P "\n";
  print P "restrictparses 0\n";
  print P "sitefile sites_ma\n";
  print P "priorfile priors_ma\n";
  print P "printsiteals 1\n";
  print P "minposterior 0.2\n";
  close P;
  return "${outMotevo}/parameters";
}

sub printLog {
  my ($outDir) = @_;
  open (L,">${outDir}/log.txt") or die "Error[startPipeline.pl]! cannot create the log ${outDir}/log.txt $!\n";
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


# sub takeBGfreq {
#   my ($a_count , $c_count , $g_count , $t_count , %bg);
#   my $tot_s = 0;
#   my $spy   = 0;
#   opendir (O,"$outDir") or die "Error[promoters_check.pl]! Cannot open $outDir\n$!\n";
#   while (my $f = readdir(O)){
#     next if (($f eq '.') || ($f eq '..'));
#     open (A ,"<${outDir}/$f") or die "Error[promoters_check.pl]! cannot open ${outDir}/$f\n";
#     while (<A>){
#       chomp $_;
#       if ($_=~/^>/){
# 	if ($_=~/$motevo_reference/){$spy = 1 ; next;}
# 	else {$spy = 0; next; }
#       }
#       #next if (($_=~/^\s*$/) or ($_=~/^>/));
#       if($spy == 1){                   #take bg frequencies of just the reference species
# 	while ($_=~/\S/g){
# 	  my $s = uc ($1);
# 	  next if ($s eq '-');
# 	  $tot_s++;
# 	  if ($s eq 'A'){$a_count++;}
# 	  elsif ($s eq 'C'){$c_count++;}
# 	  elsif ($s eq 'G'){$g_count++;}
# 	  elsif ($s eq 'T'){$t_count++;}
# 	  #else {die "Error[promoters_check.pl]! Symbol $s looks strange!\n$!\n";}
# 	}
#       }
#     }
#     close A;
#   }
#   closedir O;
#   my $a = ($a_count / $tot_s) * 100;
#   my $c = ($c_count / $tot_s) * 100;
#   my $g = ($g_count / $tot_s) * 100;
#   my $t = ($t_count / $tot_s) * 100;
#   $bg{'A'} = sprintf("%.2f",$a);
#   $bg{'C'} = sprintf("%.2f",$c);
#   $bg{'G'} = sprintf("%.2f",$g);
#   $bg{'T'} = sprintf("%.2f",$t);
#   return %bg;
# }

sub takeBGfreq {
  my (%bg);

  $bg{'A'} = 0.25;
  $bg{'C'} = 0.25;
  $bg{'G'} = 0.25;
  $bg{'T'} = 0.25;
  return %bg;
}
