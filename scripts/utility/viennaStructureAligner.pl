#!/usr/bin/perl -w

use strict;
use warnings;


if (scalar @ARGV < 2){
  die "viennaStructureAligner.pl  - given a multiple alignment and its consensus secondary structure  it aligns them
give a file with a multiple alignment and a file with its vienna format structure and it will align them.
#HINT# Here hard coded there is a oneliner to add gaps to a vienna secondary structure. The gaps come from the sequence being aligned to a MSA
";

}

my $fasta_or_aln_file  = $ARGV[0];
my $structureFile      = $ARGV[1];



#read the structure in a variable
#This works both if the structure file contains just the vienna format secondary structure (parentesis and dots)
#or if it also contains a first line and some energy value on the second line (as the standard RNAalifol output)
my ($structure , $lengthFirstLine);
open (IN,"<$structureFile") or die "cannot open the structureFile $!";
foreach my $line (<IN>){
    chomp $line;
    next if ($line =~/^\s*$/);
    $line = trim ($line);

    if (($line =~/([^\.|\(|\)])/) and (! defined $lengthFirstLine)){
      $lengthFirstLine = length ($line);
    }
    else{
      if (! defined $lengthFirstLine){
	$structure .= $line;
      }
      else{
	$structure .= substr($line, 0, $lengthFirstLine);
      }
    }
}
close IN;



my $outName = $fasta_or_aln_file;
if ($outName =~/\/([^\/]+)$/){
  $outName = $1;
}
if ($outName =~/\.aln$/){
  $outName =~s /\.aln$//;
}

print  "# STOCKHOLM 1.0\n\n";
open (IN,"<$fasta_or_aln_file") or die "cannot open the alignment File $!";
my $spy = 0;
my ($spacerSample, $spacerSize);
my ($blockSample, $blockSize);
foreach my $line (<IN>){
  chomp $line;
  next if ($line=~/CLUSTAL/);
  print  "\n" if (($line=~/^\s*$/) and ($spy == 0));

  if ($line=~/^(\S+\s+)(\S+)/){
    $spacerSample = $1;
    $blockSample  = $2;
    $spacerSize = length ($spacerSample);
    $blockSize  = length ($blockSample);
    if ($spacerSize > 14){
      print  "$line\n";
    }
    else{
      my $spacer;
      my $spaceToAdd = 14 - $spacerSize;
      my $i = 0;
      while ($i < $spaceToAdd){$i++;$spacer .= ' ';}
      print  "$spacerSample"."$spacer"."$blockSample"."\n";
    }
    $spy = 1;
    next;
  }

  if (($line=~/^\s*/) and ($spy == 1)){
    $spy  = 0;
    my $spacer;
    my $stretchStructure = substr($structure, 0, $blockSize);
    substr($structure, 0, $blockSize) = "";
    if ($spacerSize > 14){
      my $spaceToAdd = $spacerSize - 14;
      my $i = 0;
      while ($i < $spaceToAdd){$i++;$spacer .= ' ';}
      print  "#=GC SS_cons  " . "$spacer" . "$stretchStructure" . "\n";
    }
    else {
      print  "#=GC SS_cons  " ."$stretchStructure" . "\n";
    }
  }
}
print  "//";
close IN;



# remove leading or trailing whitespace characters
sub trim {
  my ($string) = @_;

  $string =~ s/^\s+//;
  $string =~ s/\s+$//;

  return $string;
}

###COMMAND LINE
#The purpose of this oneliner is to have a certain vienna format structure aligned to a fasta_aln sequence.
#This is useful in the situation in which you have a secondary structure estimated on a single sequence, without gap (like for mFOLD), and you wanna align this to the same sequence but with gaps, extracted from a Multiple Sequence Alignment
#After you run the oneliner you can append the secondary straucture generated to the MSA of interest, using viennaStructureAligner.pl

#To run the oneline you just have to edit manually the two files that are opened, the first containing the fasta_aln of the single sequence, and the other containing the vienna format secondary structure

#perl -e 'open(FA,"<pr.fasta_aln");foreach $l(<FA>){next if ($l=~/>/);chomp $l;if($l=~/(\S+)/){$allLine .= $1;}}close FA;print"$allLine\n";    $a=0;while($allLine=~/(\S)/g){$s=$1;$a++; if($s eq "-"){$gap{$a}=1;}}                                     open(SS,"<../mFold.vienna.ss");foreach $l(<SS>){next if ($l=~/[UuTtAaCcGg]/);if($l=~/^(\S+)/){$allSS=$1;}}close SS;   $a=0;while($allSS=~/(\S)/g){$s=$1;$a++;while(defined $gap{$a}){$a++;print "."}print"$s";}print"\n";'
