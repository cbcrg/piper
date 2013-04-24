#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;

my $inputDir = $ARGV[0];
my $outDir   = $ARGV[1];
my $repCov   = $ARGV[2];
if (scalar @ARGV < 3){
  die "USAGE:
        perl repeatCoverage.pl inputDirectory outputDirectory repeatCoverageThreshold
";
}


my @fastas = `ls $inputDir/*.fa`;
if ($?){die "Error [repeatCoverage.pl]! Cannot list the fasta files $!\n"};


#Take all repeat coverage
foreach my $f (@fastas){
  my %bh;
  my $name = basename($f);
  chomp $name;
  $name=~s/\.fa$//;
  system "rm ${outDir}/repeatCoverage_${name}.id"       if (-e "${outDir}/repeatCoverage_${name}.id");
  system "rm ${outDir}/${name}.bh.rep${repCov}.ex.gtf"  if (-e "${outDir}/${name}.bh.rep${repCov}.ex.gtf");
  system "rm ${outDir}/${name}.bh.rep${repCov}.fa"      if (-e "${outDir}/${name}.bh.rep${repCov}.fa");
  system "rm ${outDir}/${name}.all.rep${repCov}.ex.gtf" if (-e "${outDir}/${name}.all.rep${repCov}.ex.gtf");
  system "rm ${outDir}/${name}.all.rep${repCov}.fa"     if (-e "${outDir}/${name}.all.rep${repCov}.fa");

  open (O1, ">${outDir}/repeatCoverage_${name}.id") or die "Error [repeatCoverage.pl]! Cannot create ${outDir}/repeatCoverage_${name}.id  $! \n";
  open(F,"<$f") or die "Error [repeatCoverage.pl]! Cannot open $f $!\n";
  my ($fraction , $id , @seq);
  my $size    = 0;
  my $rep     = 0;
  my $starter = 0;
  while(my $line = <F>){
      chomp $line;
      #check
      if (($line=~/>/) && ($starter == 1)){
	  $fraction = ($rep/$size)*100;
	  print O1 "$fraction\n";
	  if ($fraction <= $repCov){
	      #print all
	      printFunc($id , $name , \@seq , "all");
	      #print best hit
	      if ($id=~/^(.+)_hit\d+$/){
		  my $tx = $1;
		  if (! defined $bh{$tx}){
		      printFunc($id , $name , \@seq , "bh");
		  }
		  $bh{$tx} = 1;
	      }
	  }
	  @seq = ();
	  $id = '';
	  $rep = 0;
	  $size= 0;
      }
      $starter = 1;
      #load
      if($line=~/>(.+)/){
	  $id = $1;
	  print O1 "$id\t";
      }
      else {
	  push (@seq,"$line\n");
	  while ($line=~/(\S)/g){
	      my $s = $1;
	      next if ($s eq " ");
	      $size++;
	      if(lc($s)eq($s)){
		  $rep++;
	      }
	  }
      }
  }

#dump
  $fraction = ($rep/$size)*100;
  print O1 "$fraction\n";
  if ($fraction <= $repCov){
      printFunc($id , $name , \@seq , "all");
      #print best hit
      if ($id=~/^(.+)_hit\d+$/){
	  my $tx = $1;
	  if (! defined $bh{$tx}){
	      printFunc($id , $name , \@seq , "bh");
	  }
	  $bh{$tx} = 1;
      }
  }

  close F;
  close O1;
}




sub printFunc {
    my ($id , $name , $refseq , $type) = @_;
    my @seq = @{$refseq};
    if ($type eq 'all'){
	system "grep -P \"$id\\S;\" $inputDir/${name}.ex.gtf >> ${outDir}/${name}.all.rep${repCov}.ex.gtf";
	open (O2,">>${outDir}/${name}.all.rep${repCov}.fa") or die "Error [repeatCoverage.pl]! Cannot create ${outDir}/${name}.all.rep${repCov}.fa $!\n";
	print O2 ">$id\n";
	foreach my $s (@seq){print O2 "$s";}
	close O2;
    }
    if ($type eq 'bh'){
	system "grep -P \"$id\\S;\" $inputDir/${name}.ex.gtf >> ${outDir}/${name}.bh.rep${repCov}.ex.gtf";
	open (O2,">>${outDir}/${name}.bh.rep${repCov}.fa") or die "Error [repeatCoverage.pl]! Cannot create ${outDir}/${name}.bh.rep${repCov}.fa $!\n";
	print O2 ">$id\n";
	while (my $s =shift(@seq)){print O2 "$s";}
	close O2;
    }
}
