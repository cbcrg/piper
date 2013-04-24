#!/usr/bin/perl -w

use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
die ("USAGE: Specify a transcript gtf file (es: cow.tx.gtf). It does not work with exon gtf files.
 This script will return clusters of transcripts (genes), i.e. transcript loci overlapping by at least 1nt") if scalar @ARGV<1;

###INFO: this script can be used downstream the PipeR mapping or any time you have a set of transcript loci, and you want to guess the genes. 
#These are defined as separate transcript clusters

my $input = $ARGV[0];
my $startR = -1; 
my @list   = `sort -k 1,1 -k 7,7 -k 4,4n $input`; 
my ($endR,$chrR,$strandR,$chr0,$strand0,$start0,$end0,$chr1,$strand1,$start1,$end1);

my $counter = 0;
foreach my $c (0..$#list){
    
    #first line
    if($startR == -1) {
	if ($list[$c]=~/(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+/){
	    $chr0   =$1;
	    $start0 =$2;
	    $end0   =$3;
	    $strand0=$4;
	}
    }
    else {
	$start0 =$startR;
	$end0   =$endR;
	$strand0=$strandR;
	$chr0   =$chrR;
    }  
	
    if (! defined $list[1+$c]){
	print"$chr0\ttranscriptCluster\tgene\t$start0\t$end0\t\.\t$strand0\t\.\tgene_id \"GENE_$counter\"; transcript_id \"GENE_$counter\";\n";
	next;
    }
    #second line
    if ($list[1+$c]=~/(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+/){
	$chr1    =$1;
	$start1  =$2;
	$end1    =$3;
	$strand1 =$4;
    } 
	
	
    #verify overlap
    if (($chr0 eq $chr1)&&($strand0 eq $strand1) && (($start1<=$end0)&&($start1>=$start0))){	    
	$startR =min($start0,$start1); 
	$endR   =max($end0,$end1); 
	$chrR   =$chr0; 
	$strandR=$strand0;  
	
	}
	else{
	    $counter++;
	    print"$chr0\ttranscriptCluster\tgene\t$start0\t$end0\t\.\t$strand0\t\.\tgene_id \"GENE_$counter\"; transcript_id \"GENE_$counter\";\n";
	    $startR=-1;
	}        
}
