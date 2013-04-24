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
my ($query_promoter_anchor , $splitGenome , $blastName , $queryName  , $strategyName  , $clusterName , $experimentName , $pipelineDirName) = options();
$queryName = $query_promoter_anchor unless ($query_promoter_anchor eq 'none');
my $databaseName       =  $pipelineDirName . '/allGenomeInfo/';
my $configName         = "${pipelineDirName}/experiments/${experimentName}/CONFIG/blastConfig";
my $clusterConfigName  = "${pipelineDirName}/experiments/${experimentName}/CONFIG/clusterConfig";
my $blastOutDir        = "${pipelineDirName}/experiments/${experimentName}/BLAST_OUT/";
my $clusterFileDir     = "${pipelineDirName}/experiments/${experimentName}/CLUSTER_FILES/";
my $total              = 0;


#CONFIG
my ($configLine , $clusterConfigLine , @clusterHeaderLines);
readConfig ();
readClusterConfig() if ($clusterName eq 'on');

#PREPARING
clearBlastOut();
clearDone();
opendir (D,"$databaseName") or die "Error[blastSearch.pl]! cannot read the $databaseName directory $!\n";
my @allSpeciesDir = readdir(D);
close D;
my (%dbSizes , %assemblyInfo);
if (($splitGenome eq 'yes') || ($strategyName eq 'abblastr') || ($strategyName eq 'wublastr')){
  my ($dbSizes_p , $assemblyInfo_p) = takeRealDbSize();
  %dbSizes      = %{$dbSizes_p};
  %assemblyInfo = %{$assemblyInfo_p};
}


#BLASTING
my %blastOutSingleDBnames;
foreach my $speciesDir (@allSpeciesDir){
    next if (($speciesDir eq '.') or ($speciesDir eq '..'));
    next unless (-d $databaseName . '/' . $speciesDir);
    print "#blasting against ${speciesDir}...\n";
    my ($blast_db);

    if (($strategyName eq 'wublastn') or ($strategyName eq 'wublastn_opt')){
	$blast_db  = $databaseName . '/' . $speciesDir . '/wublastn_db/';
    }
    if (($strategyName eq 'abblastn') or ($strategyName eq 'abblastn_opt')){
	$blast_db  = $databaseName . '/' . $speciesDir . '/abblastn_db/';
    }
    if ($strategyName eq 'wublastr'){
	$blast_db  = $databaseName . '/' . $speciesDir . '/wublastr_db/';
    }
    if ($strategyName eq 'abblastr'){
	$blast_db  = $databaseName . '/' . $speciesDir . '/abblastr_db/';
    }

    my @singleDBs = takeSingleDBcmd($blast_db);

    foreach my $singleDB (@singleDBs){
      chomp $singleDB;
      my $blast_cmd = $blastName . " ${blast_db}/$singleDB $queryName";
      $blast_cmd .= " -mformat=2";
      $blast_cmd .= " -Z=" . $dbSizes{$speciesDir} . " " if ($splitGenome eq 'yes');
      $blast_cmd .= " $configLine"  if (defined $configLine);
      my $outName;
      if ($splitGenome eq 'no') {$outName = "${blastOutDir}/${speciesDir}.mf2";}
      if ($splitGenome eq 'yes'){$outName = "${blastOutDir}/${speciesDir}.${singleDB}__tmp__.mf2";push (@{$blastOutSingleDBnames{$speciesDir}},$outName)}
      $blast_cmd .= " > $outName ";
      if ($clusterName eq 'on'){
	runOnCluster($blast_cmd , $speciesDir , $singleDB);
      }
      else{
	(system "$blast_cmd") == 0 or die "Error[blastSearch.pl]! Cannot execute:\n$blast_cmd\n$!\n";
      }
    }

}


#WAITING THAT ALL THE SCRIPT EXITED THE CLUSTER
waiting() if ($clusterName eq 'on');

#CLEAN _DONE_ FILES
clearDone();

clusterErrorCheck() if ($clusterName eq 'on');
if ($splitGenome eq 'yes'){
  concatenateBlastOutput();
}
sortBlastOut();
blastRmf2_to_standardMf2() if (($strategyName eq 'abblastr') || ($strategyName eq 'wublastr'));


###
#FUNCTIONS
sub options {
  my ($query_promoter_anchor , $splitGenome , $blastName , $queryName  , $strategyName  , $clusterName , $experimentName , $pipelineDirName);
  my $spyQuery_promoter_anchor = 1;
  my $spyBlast                 = 1;
  my $spyQuery                 = 1;
  my $spyStrategy              = 1;
  my $spyCluster               = 1;
  my $spyExperiment            = 1;
  my $spyPipelineDir           = 1;
  my $spySplitGenome           = 1;

  foreach my $field (0..$#ARGV){
    if ($ARGV[$field] eq '-queryPromoterAnchor'){
	$query_promoter_anchor = $ARGV[1+$field];
	$spyQuery_promoter_anchor = 2;
	next;
    }
    if ($spyQuery_promoter_anchor == 2){
	$spyQuery_promoter_anchor = 3;
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
    if ($ARGV[$field] eq '-cluster'){
	$clusterName = $ARGV[1+$field];
	$spyCluster = 2;
	next;
    }
    if ($spyCluster == 2){
	$spyCluster = 3;
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
    if ($ARGV[$field] eq '-query'){
	$queryName = $ARGV[1+$field];
	$spyQuery = 2;
	next;
    }
    if ($spyQuery == 2){
	$spyQuery = 3;
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

  die "Error[blastSearch.pl]! You must provide the -experiment parameter.\n"   if ($spyExperiment  != 3);
  die "Error[blastSearch.pl]! You must provide the -pipeline_dir parameter.\n" if ($spyPipelineDir != 3);
  die "Error[blastSearch.pl]! You must provide the -query parameter.\n"        if (! defined $queryName);
  if ((defined $strategyName) && ($strategyName ne "wublastn") && ($strategyName ne "abblastn") && ($strategyName ne "wublastr") && ($strategyName ne "abblastr") && ($strategyName ne "abblastn_opt") && ($strategyName ne "wublastn_opt")){
      die "Error[blastSearch.pl]! The supported blast flavors are wublastn abblastn blastr abblastn_opt wublastn_opt\n";
  }
  die "Error[blastSearch.pl]! Allowed cluster options are on|off \n"                                                if ((defined $clusterName) && ($clusterName ne 'off') && ($clusterName ne 'on'));
  die "Error[blastSearch.pl]! -splitGenome parameter can be either \'yes\' or \'no\'\n"                             if ((defined $splitGenome) && ($splitGenome ne 'yes') && ($splitGenome ne 'no'));
  die "Error[blastSearch.pl]! The -splitGenome function is meant to run on the cluster, please set -cluster on\n"   if ((defined $splitGenome) && ($splitGenome eq 'yes') && ($clusterName ne 'on'));
  $blastName             = 'wu-blastn'    if (! defined $blastName);
  $strategyName          = 'wublastn_opt' if (! defined $strategyName);
  $clusterName           = 'off'          if (! defined $clusterName);
  $splitGenome           = 'no'           if (! defined $splitGenome);
  $query_promoter_anchor = 'none'         if (! defined $query_promoter_anchor);
  return ($query_promoter_anchor , $splitGenome , $blastName , $queryName  , $strategyName  , $clusterName , $experimentName , $pipelineDirName);
}

sub acceptedVariableSpace {
  my %space = ('-queryPromoterAnchor' => 1 , '-splitGenome' => 1 , '-cluster' => 1 , '-blast' => 1  , '-query' => 1  , '-blast_strategy' => 1  , '-experiment' => 1  , '-pipeline_dir' => 1 );
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[blastSearch.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}


sub help_message {
my $helpMessage = "\nNAME
blastSearch.pl - Perform the blast search using the queries and the target genomes defined in the selected experiment\n
SYNOPSIS
blastSearch.pl -experiment -pipeline_dir -query [-cluster [on|off] -blast -blast_strategy] \n
DESCRIPTION
   * blastSearch.pl takes as input a query multi-FASTAfile the experiment name and the pipeline directory
   * it returns the output (in mformat=2) in the \"BLAST_OUT\" experiment folder .
OPTIONS
   * By default blastSearch.pl assumes that blast program name on your computer is \'wu-blastn\'.
     If this is not the case you must provide the field -blast with the appropriate path.
     If the user want to use abblast or blastr must specify with -blast the proper blast name
   * By default blastSearch.pl assumes that the blast strategy that will be used is wublastn_opt.
     You can change this by editing the -strategy parameter. The supported strategies are: wublastn|abblastn|abblastn_opt|wublastn_opt|wublastr|abblastr
   * The user can run the blast versus the individual target genomes on the cluster (if available). BlastSearch.pl will generate bash scripts with the blast commandline in the CLUSTER_FILES experiment folder. 
     The cluster will be called using the command line specified under \"#COMMAND\" in CONFIG/clusterConfig, while the bash scripts will include the header lines specified under \"#COMMAND\" in CONFIG/clusterConfig. 
     BlastSearch.pl will wait untill all the jobs return an exit code from the cluster.

TROUBLESHOOTING
   * Make sure that the blast you use in in agreement with the formatting of the database.
   * The user can edit the CONFIG/blastConfig in the experiment folder. This file is ment to set the blast parameters (i.e. e-values thresholds...) to be used. By default an e-value threshold of 0.00001 and the parameters specific of each blast strategy are used.

";
print "$helpMessage\n\n\n";
exit;
}


sub readConfig {
    open (C,"<$configName") or die "Error[blastSearch.pl]! Cannot read the blast configuration file $!\n";
    foreach my $line (<C>){
	chomp $line;
	next if ($line=~/^\s*$/);
	$configLine .= $line . " ";
    }
    if (($splitGenome eq 'yes') and (($configLine =~/-Z/) or ($configLine =~/Z=/))){die "Error[blastSearch.pl]! Please erease from the blastConfig the Z option. This is not allowed when running on splitGenome mode \n";}
    if ($configLine =~/mformat/){ die "Error[blastSearch.pl]! Please erease from the blastConfig the mformat option \n";}
    close C;
}

sub readClusterConfig {
    open (CC,"<$clusterConfigName") or die "Error[blastSearch.pl]! Cannot read $clusterConfigName \n$!\n";
    my $spyHeader = 0;
    foreach my $line (<CC>){
	chomp $line;
	next if ($line=~/^\s*$/);
	$clusterConfigLine = $line if (($line=~/##SCRIPT##/) and ($spyHeader == 0));

	if ($line=~/^#HEADER/){$spyHeader = 1;}
	push (@clusterHeaderLines , $line) if ($spyHeader == 1);
    }
    close CC;
}



sub runOnCluster {
    my ($blast_cmd , $speciesDir , $singleDB) = @_;

    #create the script
    my $scriptName = "${clusterFileDir}/${speciesDir}.${singleDB}.bl.sh";
    open (S,">$scriptName") or die "Error[blastSearch.pl]! Cannot create the cluster script for $speciesDir $!\n";
    foreach my $h (@clusterHeaderLines){
	print S "$h\n";
    }
    print S "$blast_cmd" . "\n\n";
    print S "touch ${clusterFileDir}/${speciesDir}.${singleDB}.bl.___blast_done___" . "\n";
    close S;

    #do the cmd
    my $cmd = $clusterConfigLine;
    $cmd =~s/##SCRIPT##/$scriptName/;
    (system "$cmd") == 0 or die "Error[blastSearch.pl]! cannot run on the cluster:\n$cmd\n$!\n";
    $total++;

}

sub waiting {
    my $time = 15;
    my $done = 0; #die "$total $done\n";
    while ($done < $total ){
	print STDERR "#Waiting the job to exit the cluster... " . $done . " out of " . $total . " already finished\n";
	system "sleep $time";
	#$time += 20;
	$done = `ls ${clusterFileDir}/*___blast_done___ 2> /dev/null | wc -l`;
	chomp $done;
    }
    print STDERR "#All jobs exited the cluster\n";
}

sub clusterErrorCheck {
  my (@filesHavingError , $error , $errOut);
  my @errors = ('alloc' , 'ERROR' , 'error' , 'Error' , 'Killed' , 'memory');
  my $spy_err = 0;
  opendir (CD,"$clusterFileDir") or die "Error[blastSearch.pl]! cannot open the directory $clusterFileDir\n$!\n";
  while( ( my $file = readdir(CD))){
    next if (($file eq '.') or ($file eq '..') or ($file !~/\.bl\.sh\.e\S+$/));
    my $filename = "${clusterFileDir}/$file";
    foreach $error (@errors){
      $errOut = `grep -l $error $filename`;
      if ($errOut) {push (@filesHavingError , $errOut );$errOut = ''; $spy_err = 1; last;}
      $errOut = '';
    }
    if ($spy_err == 0){
      (system "cat $filename >> ${pipelineDirName}/experiments/${experimentName}/STDERR/blast") == 0 or die "Error[blastSearch.pl]! concatenate blast error message for $filename in STDERR $!\n";
      (system "rm $filename") == 0 or die "Error[blastSearch.pl]! cannot remove $filename\n$!\n";
    }
    $spy_err = 0;
  }
  closedir CD;
  if (@filesHavingError){
    print "Error[blastSearch.pl]! The following files contain the word \"error\" or \"malloc\"  or \"Killed\" or \"memory\" while running on the cluster:\n";
    foreach my $f (@filesHavingError){
      print $f ;
    }
    die "Please rerun manually the individual shell scripts before listed\n";
  }
}


sub clearDone {
    my $check_blast_done_ =  `ls ${clusterFileDir}/*___blast_done___ 2> /dev/null | wc -l`;
    chomp $check_blast_done_;
    if ($check_blast_done_ > 0){
    	system "rm ${clusterFileDir}/*___blast_done___";
    }
}

sub clearBlastOut {
    my $check_blast_out =  `ls $blastOutDir 2> /dev/null | wc -l`;
    chomp $check_blast_out;
    if ($check_blast_out > 0){
    	system "rm ${blastOutDir}/*";
    }
}

sub sortBlastOut {
  opendir (BOF , "$blastOutDir") or die "Error[blastSearch.pl]! cannot open $blastOutDir to sort\n$!\n";
  my @allBlastOutFiles = readdir(BOF);
  closedir BOF;
  foreach my $file (@allBlastOutFiles){
    next if (($file eq '.') or ($file eq '..'));
    (system "sort -k1 -k2 -k3 -k4 -k5 -k6 -k7 -k8 -k9 -k10 -k11 -k12 -k13 -k14 -k15 -k16 -k17 ${blastOutDir}/${file} > $blastOutDir/${file}_2") == 0 or die "Error[blastSearch.pl]! cannot sort\nsort -k1 -k2 -k3 -k4 -k5 -k6 -k7 -k8 -k9 -k10 -k11 -k12 -k13 -k14 -k15 -k16 -k17 ${blastOutDir}/${file} > $blastOutDir/${file}_2\n$!\n";
    (system "mv $blastOutDir/${file}_2 ${blastOutDir}/${file}")         == 0 or die "Error[blastSearch.pl]! cannot move\nmv $blastOutDir/${file}_2 ${blastOutDir}/${file}\n$!\n";
  }

}

sub concatenateBlastOutput {
  foreach my $speciesDir (@allSpeciesDir){
    next if (($speciesDir eq '.') or ($speciesDir eq '..'));
    foreach my $file (@{$blastOutSingleDBnames{$speciesDir}}){
      (system "cat $file >> ${blastOutDir}/${speciesDir}.mf2") == 0 or die "Error[blastSearch.pl]! Cannot concatenate file $file\n$!\n" ;
      (system "rm $file") == 0 or die "Error[blastSearch.pl]! Cannot remove file $file\n$!\n" ;
    }
  }
}

sub takeRealDbSize {
  my (%dbSizes , %assemblyInfo);
  foreach my $speciesDir (@allSpeciesDir){
    next if (($speciesDir eq '.') or ($speciesDir eq '..'));
    my $size         = 0;
    my $assembly_dir = "${databaseName}/${speciesDir}/chr/";
    opendir ( DIR, $assembly_dir ) || die "Error[blastSearch.pl]! cannot open dir $assembly_dir\n";
    my @allFiles = readdir(DIR);
    closedir(DIR);
    foreach my $file (@allFiles){
      next if (($file eq '.') or ($file eq '..'));
      open (F,"<${assembly_dir}/$file") or die "Error[blastSearch.pl]! cannot read ${assembly_dir}/$file\n";
      my $name;
      my $sizeXchr = 0;
      foreach my $line (<F>){
	chomp $line;
	$line =~ s/ $//g;
	next if ($line=~/^\s*$/);
	if ($line=~/>(\S+)/){
	  $name = $1;
	  next;
	}
	$sizeXchr += length($line);
	$size += length($line);
      }
      close F;
      $assemblyInfo{$speciesDir}{$name} = $sizeXchr;
    }
    $dbSizes{$speciesDir} = $size;
  }
  return (\%dbSizes , \%assemblyInfo);
}

sub blastRmf2_to_standardMf2 {
    foreach my $speciesDir (@allSpeciesDir){
      next if (($speciesDir eq '..') or ($speciesDir eq '.'));
	open (MF2 , "<${blastOutDir}/${speciesDir}.mf2") or die "Error[blastSearch.pl]! cannot open ${blastOutDir}/${speciesDir}.mf2\n$!\n";
	open (MF2new , ">${blastOutDir}/${speciesDir}.mf2.new") or die "Error[blastSearch.pl]! cannot open ${blastOutDir}/${speciesDir}.mf2.new\n$!\n";
	while (<MF2>) {
	    my $line = $_;

	    #fix the 1 nucleotide shift to the target due to the dinucleotides
	    if ($line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$/){
	      my $q_x   = $1;
	      my $q_y   = $2;
	      my $frame = $3;
	      my $t_x   = $4;
	      my $t_y   = $5;
	      unless (($q_x == 1) or ($t_x == 1)){
		$q_x = $q_x + 1;
		$t_x = $t_x + 1;
	      }
	      $q_y = $q_y + 1;
	      $t_y = $t_y + 1;
	      $line=~s/\S+\s+\S+\s+\S+\s+\S+\s+\S+$/$q_x\t$q_y\t$frame\t$t_x\t$t_y/;
	    }

	    #revers the mapping coordinates if mapping on the reverse strand
	    if($line=~/^\S+\s+(\S+)__rev__\s+/){
	      my $chr    = $1;
	      my ($chrLen , $x , $y , $x_false , $y_false);
	      if (defined $assemblyInfo{$speciesDir}{$chr}){
		$chrLen = $assemblyInfo{$speciesDir}{$chr};
	      }
	      else {
		die "Error[blastSearch.pl]! chr $chr not defined!\n";
	      }
	      if ($line=~/(\S+)\s+(\S+)$/){
		$x_false = $1;
		$y_false = $2;
		$x       = ($chrLen - $y_false) +1;
		$y       = ($chrLen - $x_false) +1;
	      }
	      $line=~s/\S+\s+\S+$/$x\t$y/;
	    }

	    #rewrite the strands
	    if ($line=~/\+0(\s+\S+\s+\S+\s+)\+0(\s+\S+\s+\S+\s*)$/){
		my $f1 = $1;
		my $f2 = $2;
		$line=~s/\+0\s+\S+\s+\S+\s+\+0\s+\S+\s+\S+\s*$/\+1$f1\+1$f2/g;
	    }
	    if ($line=~/^(\S+\s+\S+)__rev__\s+/){
		my $i = $1;
		$line=~s/${i}__rev__/$i/;
		if ($line=~/\+1(\s+\S+\s+\S+\s*)$/){
		    my $v = $1;
		    $line=~s/\+1\s+\S+\s+\S+\s*$/-1$v/;
		}
	    }
	    print MF2new $line ;
	}
	close MF2;
	close MF2new;
    }

    foreach my $speciesDir (@allSpeciesDir){
      next if (($speciesDir eq '..') or ($speciesDir eq '.'));
	(system "mv ${blastOutDir}/${speciesDir}.mf2.new ${blastOutDir}/${speciesDir}.mf2") == 0 or die "Error[blastSearch.pl]! Cannot do\nmv ${blastOutDir}/${speciesDir}.mf2.new ${blastOutDir}/${speciesDir}.mf2\n$!\n";
    }
}

sub takeSingleDBcmd {
  my ($blast_db) = @_;
  my %dbs;
  opendir (DB,"$blast_db") or die "Error[blastSearch.pl]! Cannot open $blast_db\n$!\n";   #readdir
  while( (my $file = readdir(DB))){
    next if (($file eq '.') or ($file eq '..'));
    $file =~s/\.\S\S\S$//;
    $dbs{$file} = 1;
  }
  closedir (DB);
  my @singleDBs = keys(%dbs);
  undef %dbs;
  return @singleDBs;
}
