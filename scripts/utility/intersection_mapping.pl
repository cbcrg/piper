#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

if (scalar @ARGV < 3 ) {
die "intersection_mapping.pl
This script is used to compare the mapping done by different pipeR methods.
you need to specify just EXON gtf files.
The first argument is the ex.gtf file representing the query file, including all the considered transcripts
Then you have to specify the ex.gtf files of the other mapping methods (at least 2).
It is important that such files have different names, perhaps specifying the methods names (i.e. blastN_mouse.all.rep20.ex.gtf , blastR_mouse.all.rep20.ex.gtf)

OUTPUT: This script returns 2 files:
        1) txFoundByMoreMethods: list the query IDs and the methods that managed to return a successful mapping
        2) homologsSupportedByMoreMethods: each line reports the homolog_tag and the methods that managed to find exactly that homolog_tag.
           An homolog_tag is nothing but the sorted and concatenated list of exons (chr_start_end_strand;) of each homolog, plus the query_id at the end.
           In this analysis all the homologs detected by all queries and all methods, are first converted in uniq homolog_tag, then are pooled together,  
           then are checked the methods that posses that specific homolog_tag. Methos are separated by \";;;\".

        In both output files, if you just want the list of methods (i.e. to build a ven diagram) use cut -f 2 -d \" \" outFile
";
}

my (%all_tx_id , %all_methods , %toSort);

#read the query file to take the transcript IDs
open (F,"<$ARGV[0]") or die "error! cannot open the $ARGV[0] file\n$!\n";
while (my $l = <F>){
    if ($l =~/transcript_id \"([^\"]+)\"/){
	my $id = $1;
	$all_tx_id{$id} = 1;
    }
}
close F;



#transcripts found by each method
foreach my $i (1..$#ARGV){
  my $fileName = $ARGV[$i];
  my $method = basename ($fileName);
  $all_methods{$method}=();

  open (F,"<$fileName") or die "error! cannot open $fileName $!\n";
  while (my $line=<F>){
    my ($txId,$hitId,$tag,$chr,$start,$end,$strand);
    if ($line=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+/){
      $chr    = $1;
      $start  = $2;
      $end    = $3;
      $strand = $4;
      $tag = $chr . "_" . $start . "_" . $end . "_" . $strand;
      #print "$tag\n";
    }

    if ($line=~/transcript_id \"([^\"]+)\"/){
      $txId = $1;
      $all_methods{$method}{'txFound'}{$txId} = 1;
    }
    if ($line=~/hitName \"([^\"]+)\"/){
      $hitId = $1;
      $toSort{$method}{$txId}{$hitId}{$start} = $tag;
    }
  }
  close F;
}

#sanity check
die "Error! the exon gtf file you provide should have different names! Otherwise it is impossible distinguihing between the different mapping strategies\n" if ((scalar (@ARGV) -1) > scalar keys %all_methods );


#sort the exon tags and generate a uniqu tag representing the entire hit
my %hits;
foreach my $method (keys %toSort){
  foreach my $txId (keys %{$toSort{$method}}){
    foreach my $hitId (keys %{$toSort{$method}{$txId}} ){
      #foreach my $start (  keys %{$toSort{$method}{$txId}{$hitId}}){
      foreach my $start ( sort {$a <=> $b}  keys %{$toSort{$method}{$txId}{$hitId}}){
	#print "## $method  $txId  $hitId $start $toSort{$method}{$txId}{$hitId}{$start}\n";
	my $exonTag = $toSort{$method}{$txId}{$hitId}{$start};
	$hits{$method}{$txId}{$hitId} .=  "$exonTag;"  ;
      }
      $hits{$method}{$txId}{$hitId} .= "___$txId" ;
    }
  }
}



#check tx found by more methods
open (O1 , ">txFoundByMoreMethods") or die "cannot create txFoundByMoreMethods\n$!\n";
foreach my $tx (keys %all_tx_id) {
    print O1 "$tx ";
    my @tmp;
    foreach my $method (sort (keys %all_methods)){
	push (@tmp , $method) if (defined $all_methods{$method}{'txFound'}{$tx});
    }
    @tmp = sort (@tmp);
    foreach my $m (@tmp) {
	print O1 "$m;;;";
    }
    print O1 "\n";
}
close O1;


#take all homologous hits
my %allHomologs;
foreach my $txId (keys %all_tx_id) {
    foreach my $method (sort (keys %hits)){
	foreach my $hitId (keys %{$hits{$method}{$txId}} ){
	    my $tag = $hits{$method}{$txId}{$hitId};
	    $allHomologs{$tag} = 1;
	}
    }
}


#homologs supported by more methods
open (O2,">homologsSupportedByMoreMethods") or die "Error! cannot create homologsSupportedByMoreMethods $!\n";
foreach my $T (keys %allHomologs){
    print O2 "$T ";
    foreach my $txId (keys %all_tx_id) {
	foreach my $method (sort (keys %hits)){
	    foreach my $hitId (keys %{$hits{$method}{$txId}} ){
		my $tag = $hits{$method}{$txId}{$hitId};
		if ($T eq $tag){print O2 "$method".";;;" ;}
	    }
	}
    }
    print O2 "\n";
}
close O2;


# homologs supported by more methods
#foreach my $txId (keys %all_tx_id) {
#  my %foundSingle;
#  my %foundMultiple;
#  foreach my $method1 (keys %hits){
#     $foundSingle{$method1} = 1 if (defined $hits{$method1}{$txId});
#     foreach my $hitId1 (keys %{$hits{$method1}{$txId}} ){
#       my $hitTag1 = $hits{$method1}{$txId}{$hitId1};
#       my %seen;
#       foreach my $method2 (keys %hits){
#	 #skipp
#	 next if (defined $seen{$method1}{$method2});
#	 $seen{$method1}{$method2} = 1;
#	 $seen{$method2}{$method1} = 1;
#	 next if ($method1 eq $method2);
#	 foreach my $hitId2 (keys %{$hits{$method2}{$txId}} ){
#	   my $hitTag2 = $hits{$method2}{$txId}{$hitId2};
#	   if ($hitTag1 eq $hitTag2){
#	     $foundMultiple{$method1} = 1;
#	     $foundMultiple{$method2} = 1;
#	     delete( $foundSingle{$method1} );
#	     last;
#	   }
#	 }
#       }
#     }
#   }
# }


# foreach my $txId (keys %all_tx_id) {
#   my %tmp;
#   foreach my $method ($all_methods{$method}){
#     $tmp{$method} = ();

#   }

# }
