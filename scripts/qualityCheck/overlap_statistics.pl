#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use File::Basename;
#Author: Giovanni Bussotti

#HELP
my $help  = 0;
foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '-help') or ($ARGV[$field] eq '--help'));
}
if ($help > 0){
 help_message();
}

#GLOBAL VARIABLES
my ($tx_gtf , @allGtfMap ,  %ordered_ex , @out_realOverlappingFeatures , @out_tmpOverlapEx  , @out_tmpOverlapTx ,  @out_tmpOverlapNt , $extension_on_gtf2 , $pickUpGtf1Randomly , %gtf_annotations_intervalSizes , @block_sizes , @deltas4block , @blocks , @blockStarts );
my ($type , $repeatThreshold , $stranded , $pickUp_gtf_map_randomly , $extension_on_gtf_annotations , $gtf_annotations , $gtfToRandomize , $randomizations , $species , $shuffleBed , $shuffleBed_gap ,  $experimentName , $pipelineDirName) = options();
acceptedVariableSpace();
my $assemblyDir = "${pipelineDirName}/allGenomeInfo/${species}/chr/";
my $ovName      = basename($gtf_annotations);
my $outDir      = "${pipelineDirName}/experiments/${experimentName}/qualityCheck/overlapStatistics/${ovName}_VS_${species}_${type}_rep$repeatThreshold/";
(system "mkdir -p $outDir") == 0 or die "Error[overlap_statistics.pl]! Cannot execute mkdir -p $outDir\n";
printLog($outDir);
open (OUT , ">$outDir/${species}.ov.out") or die "Error[overlap_statistics.pl]! Cannot create $outDir/${species}.ov.out\n$!\n";
my $gtf_map              = "${pipelineDirName}/experiments/${experimentName}/results/${species}.${type}.rep${repeatThreshold}".'.ex.gtf';
my $hypergeometric_test  = "${pipelineDirName}/scripts/qualityCheck/hypergeometric_test.pl";

#TAKE ASSEMBLY INFO
my %infoAssembly       = takeAssemblyInfo($assemblyDir);
my $shuffleBed_chrSize = makeShufflebedChrSize(%infoAssembly) if (defined $shuffleBed);

#EXTEND GTF_ANNOTATIONS
my $ex_gtf_annotations = fileNameGenerator ('ex_gtf_annotations.gtf');
extend_real_gtf_annotations_file()     if ($extension_on_gtf_annotations > 0);
$ex_gtf_annotations = $gtf_annotations if ($extension_on_gtf_annotations == 0);

#READING GTF_MAP
if (defined $pickUp_gtf_map_randomly){
  open (GTF_MAP,"<$gtf_map") or die "Error[overlap_statistics.pl]! Cannot open $gtf_map \n$!\n";
  @allGtfMap = <GTF_MAP>;
  close GTF_MAP;
}

#REAL OVERLAPPING NTs
my $real_overlap_file = fileNameGenerator ('real_overlap_');
my $cmd_real_overlap  = "${pipelineDirName}/scripts/overlap $gtf_map $ex_gtf_annotations -m -1 -inter -o $real_overlap_file -st $stranded 2>/dev/null";
(system "$cmd_real_overlap") == 0 or die "Error[overlap_statistics.pl]! overlap failed with \n$cmd_real_overlap\n$!\n";
my %gtf_map_info_ex         = readingGTFexonsLOCAL($real_overlap_file);
my $coverage4gtf_map  = 'real';
my ($realOverlappingNts , $realOverlappingExons , $overlappingExonsList_ref , $realOverlappingTranscripts ) =  overlapping_nts (%gtf_map_info_ex);
$coverage4gtf_map    = 'random';
(system "rm $real_overlap_file") == 0 or die "Error[overlap_statistics.pl]! Cannot remove the temporary real overlap file $real_overlap_file\n$!\n";
foreach my $ov_exon_id (@{$overlappingExonsList_ref}){
  push (@out_realOverlappingFeatures, "$ov_exon_id");
}

#SORTING THE EXONS BY STARTING POSITIONS
foreach my $transcript_id (keys %gtf_map_info_ex){
   @{$ordered_ex{$transcript_id}} = sort {$gtf_map_info_ex{$transcript_id}->[$a]->{"start"} <=> $gtf_map_info_ex{$transcript_id}->[$b]->{"start"}} (0..$#{$gtf_map_info_ex{$transcript_id}});
}
%gtf_map_info_ex = ();

#TAKE EXON AND TRANSCRIPT INFO
my (%gtf_map_info_tx , $gtf_map_block , $numberTranscripts);
if ($gtfToRandomize eq 'gtf_map'){
  %gtf_map_info_ex  = readingGTFexonsLOCAL($gtf_map);
}
elsif ($gtfToRandomize eq 'gtf_annotations'){
  %gtf_map_info_ex  = readingGTFexonsLOCAL($ex_gtf_annotations);
}
%gtf_map_info_tx = gtfExons2Transcripts(%gtf_map_info_ex);
takeBlocks($gtf_map);
$gtf_map_block        = doBlocks();



#RANDOMIZATIONS
my $r               = 0;
my $overNt          = 0;
my $overExons       = 0;
my $overTranscripts = 0;
while ($r < $randomizations){
  my $ran = fileNameGenerator ('ran_');
  #creating a random gtf
  if (defined $pickUp_gtf_map_randomly){
    pickingUpfromGTFmap($ran);
  }
  elsif ($shuffleBed eq 'off'){
    createRandomGtf($ran);
  }
  else{
    createRandomWithShuffleBed($ran);
  }

  #running the overlap
  my $cmd_overlap;
  my $overlap_ran = fileNameGenerator('overlap_ran_');
  if ($gtfToRandomize eq 'gtf_map'){
    $cmd_overlap = "${pipelineDirName}/scripts/overlap $ran $gtf_annotations -m -1 -inter -o $overlap_ran -st $stranded 2>/dev/null";
  }
  if ($gtfToRandomize eq 'gtf_annotations'){
    $cmd_overlap = "${pipelineDirName}/scripts/overlap $gtf_map $ran  -m -1 -inter -o $overlap_ran -st $stranded 2>/dev/null";
  }
  (system "$cmd_overlap") == 0 or die "Error[overlap_statistics.pl]! overlap failed with\n $cmd_overlap\n";

  #read the overlap output
  my %tmp                = readingGTFexonsLOCAL("$overlap_ran");
  my ($tmp_overlappingNts , $tmp_overlappingExons , $useless , $tmp_overlappingTranscripts) =  overlapping_nts (%tmp);
  push ( @out_tmpOverlapTx , $tmp_overlappingTranscripts);
  push ( @out_tmpOverlapEx , $tmp_overlappingExons);
  push ( @out_tmpOverlapNt , $tmp_overlappingNts);

  #clean
  (system "rm $ran $overlap_ran") == 0 or die "Error[overlap_statistics.pl]! cannot clean\n$!\n";
  $r++;

  $overNt++          if ($tmp_overlappingNts         > $realOverlappingNts);
  $overExons++       if ($tmp_overlappingExons       > $realOverlappingExons);
  $overTranscripts++ if ($tmp_overlappingTranscripts > $realOverlappingTranscripts);
}


#OUTPUT
@out_tmpOverlapTx = sort {$a<=>$b} @out_tmpOverlapTx;
@out_tmpOverlapEx = sort {$a<=>$b} @out_tmpOverlapEx;
@out_tmpOverlapNt = sort {$a<=>$b} @out_tmpOverlapNt;
my ($nt_pValue , $exons_pValue , $transcripts_pValue ,  $average_tx ,  $std_tx , $zscore_tx , $median_tx , $average_ex ,  $std_ex , $zscore_ex , $median_ex , $average_nt ,  $std_nt , $zscore_nt , $median_nt);
if ($randomizations != 0){
  $nt_pValue          = $overNt    / $randomizations;
  $exons_pValue       = $overExons / $randomizations;
  $transcripts_pValue = $overTranscripts / $randomizations;
  ($average_tx ,  $std_tx , $zscore_tx , $median_tx , $average_ex ,  $std_ex , $zscore_ex , $median_ex , $average_nt ,  $std_nt , $zscore_nt , $median_nt) = statistics();
}
else{
  $nt_pValue = $exons_pValue = $transcripts_pValue = $average_tx =  $std_tx = $zscore_tx = $median_tx = $average_ex =  $std_ex = $zscore_ex = $median_ex = $average_nt =  $std_nt = $zscore_nt = $median_nt = 'NA';
}
printExtraOtputs();
print OUT "##PARAMETERS:RandomizationNumber\tExtension_on_gtf_annotations\n";
print OUT "\t$randomizations\t$extension_on_gtf_annotations\n";
print OUT "##TRANSCRIPTS:RealOverlapping\tHigerOverlap\tP-Value\tAverage\tMedian\tStandardDeviation\tZscore\n ";
print OUT "\t$realOverlappingTranscripts\t$overTranscripts\t$transcripts_pValue\t$average_tx\t$median_tx\t$std_tx\t$zscore_tx\n";
print OUT "##EXONS:RealOverlapping\tHigerOverlap\tP-Value\tAverage\tMedian\tStandardDeviation\tZscore\n ";
print OUT "\t$realOverlappingExons\t$overExons\t$exons_pValue\t$average_ex\t$median_ex\t$std_ex\t$zscore_ex\n";
print OUT "##NUCLEOTIDES:RealOverlapping\tHigerOverlap\tP-Value\tAverage\tMedian\tStandardDeviation\tZscore\n ";
print OUT "\t$realOverlappingNts\t$overNt\t$nt_pValue\t$average_nt\t$median_nt\t$std_nt\t$zscore_nt\n";
close OUT;
hypergeometricMode() if (defined $hypergeometric_test);
system "rm $ex_gtf_annotations" if ($extension_on_gtf_annotations > 0);
system "rm $shuffleBed_chrSize" if (defined $shuffleBed_chrSize);
system "rm $gtf_map_block"         if (defined $gtf_map_block);


















##FUNCTIONS
sub statistics {
  #TX
  my $total_tx = 0;
  foreach my $v (@out_tmpOverlapTx) {
    $total_tx += $v;
  }
  my $average_tx = $total_tx / @out_tmpOverlapTx;
  my $median_tx  = @out_tmpOverlapTx % 2        ? @out_tmpOverlapTx[(@out_tmpOverlapTx-1)/2]  
           :                  (@out_tmpOverlapTx[@out_tmpOverlapTx/2-1]+@out_tmpOverlapTx[@out_tmpOverlapTx/2])/2
           ;
  my $sqtotal_tx = 0;
  foreach my $v (@out_tmpOverlapTx) {
    $sqtotal_tx += ($average_tx-$v) ** 2;
  }
  my $std_tx = ($sqtotal_tx / @out_tmpOverlapTx) ** 0.5;

  my $zscore_tx;
  if ($std_tx != 0){$zscore_tx = ($realOverlappingTranscripts - $average_tx) / $std_tx;}
  else{$zscore_tx = 'NA';}

  #EX
  my $total_ex = 0;
  foreach my $v (@out_tmpOverlapEx) {
    $total_ex += $v;
  }
  my $average_ex = $total_ex / @out_tmpOverlapEx;
  my $median_ex  = @out_tmpOverlapEx % 2        ? @out_tmpOverlapEx[(@out_tmpOverlapEx-1)/2]  
           :                  (@out_tmpOverlapEx[@out_tmpOverlapEx/2-1]+@out_tmpOverlapEx[@out_tmpOverlapEx/2])/2
           ;
  my $sqtotal_ex = 0;
  foreach my $v (@out_tmpOverlapEx) {
    $sqtotal_ex += ($average_ex-$v) ** 2;
  }
  my $std_ex = ($sqtotal_ex / @out_tmpOverlapEx) ** 0.5;

  my $zscore_ex;
  if ($std_ex != 0){$zscore_ex = ($realOverlappingExons - $average_ex) / $std_ex;}
  else{$zscore_ex = 'NA';}
  #print "Min: $data[0]   Max: $data[-1]   Total: $total   count: "       . @data . "  Average: $average\n";
  #print "Median: $median   $sqtotal Standard deviation: $std\n";
  #print "Average: $average Standard deviation: $std Z-Score: $zscore";

  #NT
  my $total_nt = 0;
  foreach my $v (@out_tmpOverlapNt) {
    $total_nt += $v;
  }
  my $average_nt = $total_nt / @out_tmpOverlapNt;
  my $median_nt  = @out_tmpOverlapNt % 2        ? @out_tmpOverlapNt[(@out_tmpOverlapNt-1)/2]  
           :                  (@out_tmpOverlapNt[@out_tmpOverlapNt/2-1]+@out_tmpOverlapNt[@out_tmpOverlapNt/2])/2
           ;
  my $sqtotal_nt = 0;
  foreach my $v (@out_tmpOverlapNt) {
    $sqtotal_nt += ($average_nt-$v) ** 2;
  }
  my $std_nt = ($sqtotal_nt / @out_tmpOverlapNt) ** 0.5;

  my $zscore_nt;
  if ($std_nt != 0){$zscore_nt = ($realOverlappingNts - $average_nt) / $std_nt;}
  else{$zscore_nt = 'NA';}


  return ($average_tx ,  $std_tx , $zscore_tx , $median_tx , $average_ex ,  $std_ex , $zscore_ex , $median_ex , $average_nt ,  $std_nt , $zscore_nt , $median_nt);
}

sub random_chr {
  my $i = 0;
  my %chrNumbers;
  foreach my $chr (keys %infoAssembly){
    $chrNumbers{$i} = $chr;
    $i++;
  }
  my $range = scalar (keys %infoAssembly);
  my ($rNumGood , $rNum);
  while (! defined $rNumGood){
    $rNum = rand($range);
    $rNum=~s/\..*//;
    if ($rNum != $range){$rNumGood = 1;}
  }
  my $randomChr =  $chrNumbers{$rNum};
  return $randomChr;
}
sub random_start {
  my ($current_chr , $tr_size) = @_;
  my $range = $infoAssembly{"$current_chr"} - $tr_size;
  my $rNum = rand($range);
  $rNum=~s/\..*//;
  return ($rNum);
}
sub ovCheck {
  my ($r_chr , $r_start , $r_end , %mappedFeature) = @_;
  my $check = 1;
  foreach my $feature (0..$#{$mappedFeature{$r_chr}}){
    my $m_start = $mappedFeature{$r_chr}[$feature]{'start'};
    my $m_end   = $mappedFeature{$r_chr}[$feature]{'end'};
  if ((($m_start <= $r_end) && ($m_start >= $r_start)) || (($m_end >= $r_start) && ($m_end <= $r_end)) || (($m_start <= $r_start) && ($m_end >= $r_end)) ||   (($m_start >= $r_start) && ($m_end <= $r_end))){
      $check = 0;
      return $check;
    }
  }
  return $check;
}
sub pickingUpfromGTFmap {
  my ($ran) = @_;
  my %selection = randSubset ($pickUp_gtf_map_randomly , @allGtfMap);
  open (R,">$ran") or die "Error[overlap_statistics.pl]! cannot create the random picking up $!\n";
  foreach my $line (keys %selection){
    print R "$line";
  }
  close R;
}
sub randSubset {
  my ($subset , @allNameList) = @_;
  my $i = 0;
  my $size = scalar(@allNameList);
  my (%recall , @extractedNumbers , %outList);
  #sanity check
  if ($size < $subset){
    die "Error[overlap_statistics.pl]! the subset $subset specified is bigger than the list size $size. Impossible to extract more elements than the ones submitted! Aborting\n";
  }
  #to iterate the loop "$subset" times
  for ($i=0;$i<$subset;$i++){
    #generate a random number
    my $rNum = rand($size);
    $rNum = sprintf("%.0f", $rNum);

    #skipp already extracted numbers and the last number since the foreach below is from 0 to $size-1
    if (($recall{$rNum}) || ($rNum == $size)){
      $i--;
      next;
    }
    $recall{$rNum} = 1;
    push(@extractedNumbers , $rNum);
  }
  #generating the subset
  foreach my $index (@extractedNumbers){
    my $element = $allNameList[$index];
    $outList{$element} = 1;
  }
  return(%outList)
}

sub printExtraOtputs {
  my $distributionOfOverlappingTx = fileNameGenerator("distributionOfOverlappingTranscripts__${extension_on_gtf_annotations}__${species}__");
  open (D,">${outDir}/$distributionOfOverlappingTx") or die "Error[overlap_statistics.pl]! cannot create ${outDir}/$distributionOfOverlappingTx $!\n";
  foreach my $val (@out_tmpOverlapTx){
    print D "$val\n";
  }
  close D;
  my $distributionOfOverlappingEx = fileNameGenerator("distributionOfOverlappingExons__${extension_on_gtf_annotations}__${species}__");
  open (D,">${outDir}/$distributionOfOverlappingEx") or die "Error[overlap_statistics.pl]! cannot create ${outDir}/$distributionOfOverlappingEx $!\n";
  foreach my $val (@out_tmpOverlapEx){
    print D "$val\n";
  }
  close D;
  my $distributionOfOverlappingNt = fileNameGenerator("distributionOfOverlappingNucleotides__${extension_on_gtf_annotations}__${species}__");
  open (D,">${outDir}/$distributionOfOverlappingNt") or die "Error[overlap_statistics.pl]! cannot create ${outDir}/$distributionOfOverlappingNt $!\n";
  foreach my $val (@out_tmpOverlapNt){
    print D "$val\n";
  }
  close D;
  my $realOverlappingIds = fileNameGenerator("realOverlappingIds__${extension_on_gtf_annotations}__${species}__");
  open (O,">${outDir}/$realOverlappingIds") or die "Error[overlap_statistics.pl]! cannot create the ${outDir}/$realOverlappingIds file $!\n";
  foreach my $id (@out_realOverlappingFeatures){
    print O "$id\n";
  }
  close O;
}



#-pickUp_gtf_map_randomly' => 1 , '-extension_on_gtf_annotations' => 1 , '-gtf_annotations' => 1 ,'-gtfToRandomize' => 1 , '-randomizations' => 1 , '-experiment' => 1 , '-pipeline_dir' => 1 , '-species' => 1 , '-shuffleBed' => 1,  '-shuffleBed_gap'
#($pickUp_gtf_map_randomly , $extension_on_gtf_annotations , $gtf_annotations , $gtfToRandomize , $randomizations , $species , $shuffleBed , $shuffleBed_gap ,  $experimentName , $pipelineDirName)
sub help_message {
my $helpMessage = "\nNAME
overlap_statistics.pl - Given an annotations gff(or gtf) file (we call gtf_annotations), this script measure the overlap between the ex.gtf file of a certain species (human,cow..), a certain repeat coverage (20,80..) and mapping type (bh, all). Such file is called gtf_map. Then the input gtf_annotation (cuflink transcripts from RNAseq, chromatine markups annotations). The overlap is measured in terms of #overlapping transcripts,  #overlapping exons and #ovarlapping nucleotides. Then two kinds of statistical significance test are run:
1- Hypergeometric test considering the number of overlapping nucleotides
2- Randomization of one of the two annotation file, creating a distribution, and measuring the P-Value and the Z-Score of the real overlap (transcript exons and nucleotides) with respect to the overlap distributions

SYNOPSIS
perl overlap_statisticalRelevance.pl -gtf_annotations -experiment -pipeline_dir -species  [-randomizations -gtfToRandomize -extension_on_gtf_annotations -shuffleBed pickUp_gtf_map_randomly -shuffleBed_gap -stranded -help]

DESCRIPTION
  * give a species.ex.gtf file in the results folder and another gtf file you wanna compare with. It measures the level of overlap among the two. By default the strand of the overlap is not taken into account, just the positions are
  * It is mandatory to specify the repeat coverage of your mapping in the results folder (i.e 20 or 100)
  * It is mandatory to specify the type of your mapping in the results folder, either best hit (bh) or all.
  * If you wanna let the overlap be just on the same strand set -stranded 1
  * The script assign a significance to the gtf_map by measuring the overlap with some gtf_annotations, and assessing how much this overlap is different from the one we could get with a random gtf_map.
  * The randomization is done as in read_statistics.pl, where every time two exons belonging to different transcript overlap, the two transcript will be joined into one block. Then the blocks so defined are randomly projected on the genome. Once projected, for each transcript in each block the exons annotations are derived. Then the overlap with the gtf_annotations is measured again,    * Moreover this script considers the homologs returned by pipeR, where each different hitName is a different transcript homolog. This important because the script is going to work both when the pipeR returns for each query a single best hit (a single homolog), and when pipeR it is run exahustively, i.e. returning a one2many mapping. In the case pipeR generated multiple homologs for one query, each homolog will be considered a different transcript, with a different id. This script will tell you how much is the real overlap with provided gtf_annotations, with respect to random projections of annotations.  

  * The default projection tool is shuffleBed. However there is also a built_in projection function perfetly working. 
  * OPTION: You can decide to randomize either the gtf_annotations either the gtf_map file [Default, recommended]
  * !!!IMPORTANT!!!! The gtf_annotations can be anything, exons, enhancers, chromatine markups, gene annotations..The important is that each line contains the transcript_id field with some identifier. If you are not interested in randomizing these gtf_annotations, it really doesn't matter. But if you wanna randomize the gtf_annotations, please bear in mind that these will be reconstructed in transcript and blocks, projected, splitted again in exons. It is therefore important that the transcript_id field of each annotation is correctly assigned. On the other hand, if the annotations are indipendent one another (i.e. a set of enhancers) you should give to the gtf_annotations a transctipt_id unique for each annotation line, so that they get projected together just in the case they are overlapping.
  * OPTION: It is possible to set an extension on both side of each gtf_annotations feature to see how the statistic changes
  * The randomization are created accordingly with the chromosome assembly (-assemblyDir). The assembly dir must contain just chromosome files, one for each chromosome. (allGenomeInfo/species/chr)
  * REMARK1: The built_in projection function creates not overlapping projections (this is in fact unnecessary since the projections contiune untill the coverge is met). This function it is perfectly working, and there is an easy room for improvement. If needed it could be possible to provide a list of genomic gap were the function would not map the transcript. This is an easy implementation, copying the feature currently available for shuffleBed.
  * REMARK2: You can use -shuffleBed in order to produce randomization. This is recommendable as it can take into account regions where it is not allowed to produce randomizations, like the heterocromatine regions. If run in this way, you have to specify an extra file with te regions to exclude (via -shuffleBed_gap). These gap file are normally available for varius species via UCSC tables.
  * OPTION: -pickUp_gtf_map_randoml will randomly extract features from gtf_map (i.e. it uses real features but randomly extracted). then the rest is identical (overlap, extension...)
 * REMARK3: A bash loop is available hard copied in this script
 * The script assesses if the real overlap is significant with respect to the overlap that you find with randomly generated files.
 * The script returns all the IDs of the real overlap together with a line considering the nt coverage and the number of overlapping exons
 * OPTION: An hypergeometric test is also run to get and get the statistical relavance in this way
 * TRICK1: you can control/visualize the overlap using the display custumer track tool of UCSC http://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html#CustomTracks
 * an example of the the format UCSC accepts is hard coded to this script
 * WARNING_1!! do not run this script on multiple times on the same computer! The either the current version of overlap (by Sarah) or most probably shuffle_bed creates and delete temporary files named according with the date and the seconds of execution. If you run multiplethis script multiple time on the same machine you will have some conflicts!
 * WARNING_2!! The FPKM is thought to estimate the number of reads overlapping each transcripts. However bear in ming that mapped reads .gtf file can be enormous (i.e. 34 Giga) In this situation overlap takes forever and need an exagerated amount of memory. Use the bamtools for this purpose.

ADVICE
It makes a lot of sense and it is recommended using the default randomization of the gtf_map. You create some homology based mapping, you measure the overlap, and you generate some randomization of your mapping to see how your original mapping is far from random.
It makes a lot of sense to use an exon.gtf as gtf_annotations (i.e. generated by cufflinks after RNAseq). This is because the region of the genome that should be considered positive are the exons, not the introns. Therefore do not use transcript as gtf_annotations, but just exons.gtf.
This script will tell you the number of your homology based mapped transcript exons and nucleotides that overlap the cufflink exons. An example of command line is:
./overlap_statistics.pl -shuffleBed_gap ../../cowReadMappingData_fromDarek_UMD3_1_ens65/bosTau6_UMD_3.1_fromUCSC.gap_contig_unknown.bed -species cow -experiment exp_1 -pipeline_dir .. -gtf_annotations ../../cowReadMappingData_fromDarek_UMD3_1_ens65/mergedBt_CRGplusUS_grape_DAREK_COW_RNAseq_grapePlusCufflinks.ex.gtf

OPTIONS
  *-pipeline_dir = <path>
  exon gtf filepositin where the pipeline is installed

  *-experiment = <name>
  exon gtf filepositin where the pipeline is installedname of the experiment to run (i.e. exp_1)

  *-gtf_annotations = <file name>
  comparison gtf file. For instance an enhancer list

  *-species = <name>
  one of the species name found in the pipelineDir/allGenomeInfo   (human, for instance)

  *-repeatThreshold = <int>
  choose as gtf_map the repeat threshold from the results folder. There are always 100 80 and 80. If the user run pipeR with another repeatThreshold, he can use it

  *-type = <bh|all>
  choose as gtf_map the best hit (bh) or the all file.

  *-stranded = <0|1>
  OPTION: set it to 1 if you want to consider just the overlap happening on the same strand

  *-randomizations = <integer>
  OPTION: Number of randomization. Default 100.

  *-gtfToRandomize = <gtf_map|gtf_annotations>
  OPTION: choose the gtf to randomize (default gtf_map)

  *-extension_on_gtf_annotations = <integer>
  OPTION: Default[0]. Extend the gtf_annotations feature by a window on both sides. This extension is operated both on the real intervals and on the randomizations

  *-shuffleBed = <shufflebed program name>
  OPTION: Flag to use shuffleBed software to produce the randomization. the Default is shuffleBed. Set this parameter to 'off' if you wanna use the built_in randomization function instead.

  *-shuffleBed_gap = <file name>
  OPTION: This is mandatory if -shuffleBed. You need to specify here a .bed file indicating the regions to exclude from the shuffling. You can get this gap file for a certain assembly by UCSC tables (and selecting just chr_start_end columns)

  *-pickUp_gtf_map_randomly = <integer>
  OPTION: If you select this option the randomizations will be done by picking up a certain number of features from gtf_map (that is real features are extracted randomly)\n
";
print "$helpMessage\n\n\n";
exit;
}
sub options {
  my ($type , $repeatThreshold , $stranded , $pickUp_gtf_map_randomly , $extension_on_gtf_annotations , $gtf_annotations , $gtfToRandomize , $randomizations , $species , $shuffleBed , $shuffleBed_gap ,  $experimentName , $pipelineDirName);
  my $spyStranded                     = 1;
  my $spyPickUp_gtf_map_randomly      = 1;
  my $spyExtension_on_gtf_annotations = 1;
  my $spyGtf_annotations              = 1;
  my $spyGtfToRandomize               = 1;
  my $spyRandomizations               = 1;
  my $spySpecies                      = 1;
  my $spyShuffleBed                   = 1;
  my $spyShuffleBed_gap               = 1;
  my $spyExperiment                   = 1;
  my $spyPipelineDir                  = 1;
  my $spyRepeatThreshold              = 1;
  my $spyType                         = 1;

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
    if ($ARGV[$field] eq '-stranded'){
	$stranded = $ARGV[1+$field];
	$spyStranded = 2;
	next;
    }
    if ($spyStranded == 2){
	$spyStranded = 3;
	next;
    }
    if ($ARGV[$field] eq '-pickUp_gtf_map_randomly'){
	$pickUp_gtf_map_randomly = $ARGV[1+$field];
	$spyPickUp_gtf_map_randomly = 2;
	next;
    }
    if ($spyPickUp_gtf_map_randomly == 2){
	$spyPickUp_gtf_map_randomly = 3;
	next;
    }
    if ($ARGV[$field] eq '-extension_on_gtf_annotations'){
	$extension_on_gtf_annotations = $ARGV[1+$field];
	$spyExtension_on_gtf_annotations = 2;
	next;
    }
    if ($spyExtension_on_gtf_annotations == 2){
	$spyExtension_on_gtf_annotations = 3;
	next;
    }
    if ($ARGV[$field] eq '-gtf_annotations'){
	$gtf_annotations = $ARGV[1+$field];
	$spyGtf_annotations = 2;
	next;
    }
    if ($spyGtf_annotations == 2){
	$spyGtf_annotations = 3;
	next;
    }
    if ($ARGV[$field] eq '-gtfToRandomize'){
	$gtfToRandomize = $ARGV[1+$field];
	$spyGtfToRandomize = 2;
	next;
    }
    if ($spyGtfToRandomize == 2){
	$spyGtfToRandomize = 3;
	next;
    }
    if ($ARGV[$field] eq '-randomizations'){
	$randomizations = $ARGV[1+$field];
	$spyRandomizations = 2;
	next;
    }
    if ($spyRandomizations == 2){
	$spyRandomizations = 3;
	next;
    }
    if ($ARGV[$field] eq '-shuffleBed_gap'){
	$shuffleBed_gap = $ARGV[1+$field];
	$spyShuffleBed_gap = 2;
	next;
    }
    if ($spyShuffleBed_gap == 2){
	$spyShuffleBed_gap = 3;
	next;
    }
    if ($ARGV[$field] eq '-shuffleBed'){
	$shuffleBed = $ARGV[1+$field];
	$spyShuffleBed = 2;
	next;
    }
    if ($spyShuffleBed == 2){
	$spyShuffleBed = 3;
	next;
    }
    if ($ARGV[$field] eq '-species'){
	$species = $ARGV[1+$field];
	$spySpecies = 2;
	next;
    }
    if ($spySpecies == 2){
	$spySpecies = 3;
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
  die "Error[overlap_statistics.pl]! Please provide the -experiment parameter\n"                             if (! defined $experimentName);
  die "Error[overlap_statistics.pl]! Please provide the -pipeline_dir parameter\n"                           if (! defined $pipelineDirName);
  die "Error[overlap_statistics.pl]! Please provide the -gtf_annotations parameter\n"                        if (! defined $gtf_annotations);
  die "Error[overlap_statistics.pl]! Please provide the -species parameter\n"                                if (! defined $species);
  die "Error[overlap_statistics.pl]! Please provide the -repeatThreshold parameter (e.g. 20, 80, 100..)\n"   if (! defined $repeatThreshold);
  die "Error[overlap_statistics.pl]! Please provide the -type parameter with \"bh\" (best hit) or \"all\"\n" if (! defined $type);

  #die "Error[overlap_statistics.pl]! it is not supported using shuffleBed on the gtf_map file" if ((defined $gtfToRandomize) && ($gtfToRandomize eq 'gtf_map') && (defined $shuffleBed));

  #defaults
  $randomizations    = 100                     unless (defined $randomizations) ;
  $gtfToRandomize    = 'gtf_map'               unless (defined $gtfToRandomize) ;
  die "Error[overlap_statistics.pl]! if pickUp_gtf_map_randomly the gtf to randomize is always gtf_map\n" if     ((defined $pickUp_gtf_map_randomly)&&($gtfToRandomize eq 'gtf_annotations'));
  $extension_on_gtf_annotations = 0            unless (defined $extension_on_gtf_annotations) ;
  $stranded                     = 0            unless (defined $stranded);
  $shuffleBed                   = 'shuffleBed' unless (defined $shuffleBed);

  return ($type , $repeatThreshold , $stranded , $pickUp_gtf_map_randomly , $extension_on_gtf_annotations , $gtf_annotations , $gtfToRandomize , $randomizations , $species , $shuffleBed , $shuffleBed_gap ,  $experimentName , $pipelineDirName);
}
sub acceptedVariableSpace {
  my %space = ('-type' => 1 , '-repeatThreshold' => 1 , '-stranded' => 1 , '-pickUp_gtf_map_randomly' => 1 , '-extension_on_gtf_annotations' => 1 , '-gtf_annotations' => 1 ,'-gtfToRandomize' => 1 , '-randomizations' => 1 , '-experiment' => 1 , '-pipeline_dir' => 1 , '-species' => 1 , '-shuffleBed' => 1,  '-shuffleBed_gap' => 1 );
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[overlap_statistics.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}
sub takeAssemblyInfo {
  my ($assembly_dir) = @_;
  opendir ( DIR, $assembly_dir ) || die "Error[overlap_statistics.pl]! Cannot open dir\n $assembly_dir\n";
  my %info;
  while (my $file = readdir(DIR)){
    next if (($file eq '.') or ($file eq '..'));
    open (F,"<${assembly_dir}/$file") or die "Error[overlap_statistics.pl]! cannot read ${assembly_dir}/$file\n";
    my $name;
    my $totLen = 0;
    while (my $line = <F>){
      chomp $line;
      next if ($line=~/^\s*$/);
      $line =~ s/ //g;

      if ($line=~/>(\S+)/){
	$name = $1;
	next;
      }
      $totLen += length($line);
    }
    close F;
    $info{$name} = $totLen;
  }
  closedir(DIR);
  return %info;
}
sub fileNameGenerator{
  my ($nameRoot) = @_;
  my $tmp_name_counter = 0;
  my $tmp_name;
  while (!$tmp_name || -f $tmp_name) {
    $tmp_name_counter++;
    $tmp_name = "$nameRoot".".$tmp_name_counter";
  }
  return $tmp_name;
}
sub extend_real_gtf_annotations_file {
  open (EXT,">$ex_gtf_annotations") or die "Error[overlap_statistics.pl]!cannot create the $ex_gtf_annotations\n$!\n";
  open (G2,"<$gtf_annotations") or die "Error[overlap_statistics.pl]! cannot open the gtf_annotations file $!\n";
  while (my $line = <G2>){
    if ($line=~/^(\S+)(\t\S+\t\S+\t)(\S+)\t(\S+)(\t.+)/){
      my $chr           = $1;
      my $block1        = $2;
      my $extendedStart = (min ($3 , $4)) - ($extension_on_gtf_annotations);
      my $extendedEnd   = (max ($4 , $3)) + ($extension_on_gtf_annotations);
      my $block2        = $5;
      $extendedStart    = 1 if ($extendedStart < 0);
      $extendedEnd      = $infoAssembly{$chr} if ($extendedEnd > $infoAssembly{$chr}) ;
      print EXT "$chr$block1$extendedStart\t$extendedEnd$block2\n";
    }
  }
  close G2;
  close EXT;
}
sub readingGTFexonsLOCAL {
  my ($GTFfile) = @_;
  open (IN,"<$GTFfile") or die "Error[overlap_statistics.pl]! cannot open the GTFfile $!";
  my %infoGtf;
  foreach my $line (<IN>){
	chomp $line;
	my ($chr , $group , $frame , $strand , $score , $end , $start , $feature , $source , $gene_id , $transcript_id);
	if ($line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.*)$/){
	    $chr     = $1;
	    $source  = $2;
	    $feature = $3;
	    $start   = min ($4, $5);
	    $end     = max ($4, $5);
	    $score   = $6;
	    $strand  = $7;
	    $frame   = $8;
	    $group   = $9;
	    if ($group =~/gene_id \"([^\"]+)\"/){
		$gene_id = $1;
	    }
	    else {die "Error[overlap_statistics.pl]! gtf file in wrong format. Impossible to find the gene_id in line:\n$line\n when it is mandatory for this format\n";}
	    if ($group =~/hitName \"([^\"]+)\"/){
		$transcript_id = $1;
	    }
	    else {die "Error[overlap_statistics.pl]! gtf file in wrong format. Impossible to find the transcript_id in line:\n$line\n when it is mandatory for this format\n";}

	   # next unless ($feature eq 'exon');

	    my %hash = (
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
	else {die "Error[overlap_statistics.pl]! cannot read the line $line\n";}
      }
  close IN;
  return %infoGtf;
}
sub overlapping_nts {
  my (%info) = @_;
  my (%hash , @exonIdList);
  my $ntOverlap4set         = 0;
  my $exonOverlap4set       = 0;
  my $transcriptOverlap4set = 0;
  foreach my $transcript_id (keys %info){
    my $overlap4transcript          = 0;
    my $numberOfOverlappingFeatures = 0;
    my $numberOfOverlappingExons    = 0;
    foreach my $index (0..$#{$info{$transcript_id}}){
      my $exonGroup = $info{$transcript_id}->[$index]->{'group'};
      if ($exonGroup=~/nb_ov_feat2:\s+(\d+)\s+list_feat2:\s+(.*)/){
	my $numberOfFeatures = $1;
	my $rest             = $2;
	if ($numberOfFeatures > 0){
	  while ($rest){
	    if ($rest=~/^.+_([^_]+)_([^_]+)_\.,(.*)$/){
	      my $start = min($1,$2);
	      my $end   = max($1,$2);
	      my $currentOverlappingNucleotides = abs($end - $start) + 1;
	      $rest  = $3;
	      $overlap4transcript += $currentOverlappingNucleotides;
	    }
	    else{
	      $rest='';
	    }
	  }
	  $numberOfOverlappingFeatures += $numberOfFeatures;
	  $numberOfOverlappingExons++;
          push(@exonIdList,$transcript_id);
	}
      }
    }
    $hash{$transcript_id}{'overlappingNucleotides'}      = $overlap4transcript;                        #foreach gtf_map transcripts, this is the sum of gtf_annotation NTs that overlap
    $hash{$transcript_id}{'numberOfOverlappingFeatures'} = $numberOfOverlappingFeatures;               #foreach gtf_map transcripts, this is the sum of gtf_annotation that overlap
    $hash{$transcript_id}{'numberOfOverlappingExons'}    = $numberOfOverlappingExons;                  #exons of gtf_map transcripts that are overlapped by at least 1 gtf_annotation
    $ntOverlap4set   += $overlap4transcript;
    $exonOverlap4set += $numberOfOverlappingExons;
    $transcriptOverlap4set++ if ($numberOfOverlappingExons > 0);
  }

  #MAKE COVERAGE FILE
  if ($coverage4gtf_map eq 'real'){
    my %mature_tx_length;
    my $numberOfFeatures = `wc -l $gtf_annotations | cut -f1 -d " " | tr -d "\n"`; if ($?){die "Error[overlap_statistics.pl]! Cannot count the number if features in $gtf_annotations\n";}
    my $millionOfFeatures = $numberOfFeatures / 1000000;

    my $coverage4tx = fileNameGenerator("coverage4${species}_tx");
    open (C,">${outDir}/$coverage4tx") or die "Error[overlap_statistics.pl]! Cannot create ${outDir}/$coverage4tx\n$!\n";
    print C "###${species}_tx_ID\t#overlapping_nt\t#overlapping_events\t#exons_having_at_least_1overlapping_feature\tFeatures_Per_Kilobase_of_transcript_model_per_Million_mapped_Features(FPKM)\n";
    foreach my $tx_id (keys %hash){
      $mature_tx_length{$tx_id} = 0;
      foreach my $index (0..$#{$info{$tx_id}}){
	$mature_tx_length{$tx_id} += abs($info{$tx_id}->[$index]->{'end'} - $info{$tx_id}->[$index]->{'start'}) + 1 ;
      }
      my $FPKM = sprintf ("%.0f" , $hash{$tx_id}{'numberOfOverlappingFeatures'} / ($millionOfFeatures * ($mature_tx_length{$tx_id} / 1000))) ;
      print C "$tx_id" . "\t" . $hash{$tx_id}{'overlappingNucleotides'} . "\t" . $hash{$tx_id}{'numberOfOverlappingFeatures'} . "\t" . $hash{$tx_id}{'numberOfOverlappingExons'} . "\t" . $FPKM . "\n";
    }
    close C;
  }

  return ( $ntOverlap4set , $exonOverlap4set , \@exonIdList , $transcriptOverlap4set );
}
sub makeShufflebedChrSize {
  my (%infoAssembly) = @_;
  my $shuffleBed_chrSize = fileNameGenerator('shuffleBed_chrSize');
  open (S,">$shuffleBed_chrSize") or die "Error[overlap_statistics.pl]! Cannot create $shuffleBed_chrSize $!\n";
  foreach my $chr (keys %infoAssembly){
    print S "$chr\t" .  $infoAssembly{$chr} . "\n";
  }
  close S;
  return $shuffleBed_chrSize;
}

sub verifyOverlap {
  my ($traSta ,  $traEnd , $geneSta , $geneEnd ) = @_;
  my $overlap = 0;
  if ((($traSta <= $geneEnd) && ($traSta >= $geneSta)) || (($traEnd >= $geneSta) && ($traEnd <= $geneEnd)) || (($traSta <= $geneSta) && ($traEnd >= $geneEnd)) ||   (($traSta >= $geneSta) && ($traEnd <= $geneEnd))){
    $overlap = 1;
  }
  return $overlap;
}
sub takeCoverage {
  my ($current_gtf) = @_;
  my $cmd_sort = 'sort -k 1,1 -k 7,7 -k 4,4n  ' . "$current_gtf";
  my @sorted_gtf = `$cmd_sort`; if ($?) {die "Error[overlap_statistics.pl]! cannot run $cmd_sort\n";}

  my $cov = 0;
  my $l   = shift( @sorted_gtf);
  my ($oldChr , $oldStart , $oldEnd ,$oldStrand);
  if ($l=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+/){
       $oldChr    = $1;
       $oldStart  = $2;
       $oldEnd    = $3;
       $oldStrand = $4;
  }

  while (my $l = shift(@sorted_gtf)){
      if ($l=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+/){
	  my $chr    = $1;
	  my $start  = $2;
	  my $end    = $3;
	  my $strand = $4;

	  if (($chr eq $oldChr) && ($strand eq $oldStrand)){
	      my $ov = verifyOverlap($oldStart , $oldEnd , $start , $end);
	      if ($ov == 1) {
		  $oldStart = min ($oldStart , $start);
		  $oldEnd   = max ($oldEnd , $end);
		  next;
	      } 
	      else{
		  $cov += abs($oldEnd - $oldStart) + 1;
		  $oldStart = $start; $oldEnd = $end;
		  next;
	      }
	  }
	  else{
	      $cov += abs($oldEnd - $oldStart) + 1;
	      $oldChr = $chr; $oldStart = $start; $oldEnd = $end; $oldStrand = $strand;
	      next;
	  }
      }
  }
  $cov += abs($oldEnd - $oldStart) + 1;
  return $cov;
}


sub transcriptGtf2exonGtf {
    my ($in_tx_file , $gtf_map_info_ex_ref) = @_;
    my %ex_hash = %{$gtf_map_info_ex_ref};
    my %out;
    open (IN,"<$in_tx_file") or die "Error[overlap_statistics.pl]! cannot read $in_tx_file\n$!\n";
    while (my $line = <IN>){
	my ($chr , $start , $strand , $tx_id);
	if ($line =~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+/){
	    $chr    = $1;
	    $start  = $2;
	    $strand = $3;
	}
	if ($line =~/hitName "([^\"]+)"/ ){
	    $tx_id = $1;
	}
	my $delta = -5;
	my ($newStart , $newEnd);
	foreach my $ex (sort {$a->{"start"} <=> $b->{"start"}}  @{$ex_hash{$tx_id}}){
	    $delta = abs($ex->{"start"} - $start) if ($delta == -5);
	    if ($start > $ex->{"start"}){
		$newStart = $ex->{"start"} + $delta;
		$newEnd   = $ex->{"end"}   + $delta;
	    }
	    if ($start <= $ex->{"start"}){
		$newStart = $ex->{"start"} - $delta;
		$newEnd   = $ex->{"end"}   - $delta;
	    }
	    my $newExon = "$chr\t" . $ex->{"source"} . "\t" . $ex->{"feature"} . "\t$newStart\t$newEnd\t" . $ex->{"score"} . "\t$strand\t" . $ex->{"frame"} . "\t" . $ex->{"group"};
	    push (@{$out{$tx_id}},$newExon);
	}
    }
    close IN;
    return %out;
}
sub gtfExons2Transcripts{
  my (%infoGtf) = @_;
  my %infoGtfTranscript;
  foreach my $rna (keys %infoGtf){
    my $minExonStart = 1000000000000000000000000000000000000000000000000000000000000000000000000000;
    my $maxExonEnd   = 0;
    my ($chr , $group , $frame , $strand , $score , $feature , $source , $gene_id);
    foreach my $index (0..$#{$infoGtf{$rna}}){
      $minExonStart = min ($minExonStart, $infoGtf{$rna}->[$index]->{"start"});
      $maxExonEnd   = max ($maxExonEnd, $infoGtf{$rna}->[$index]->{"end"});

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
 sub hypergeometricMode {
   #take the population size
   my $allSize = 0;
   foreach my $current_chr (keys %infoAssembly){
     $allSize += $infoAssembly{"$current_chr"};
   }

   #take the size of the two compared subgroups (in nts coverage)
   my $count1 = 0;
   foreach my $transcript_id (keys %gtf_map_info_ex){
     foreach my $index (0..$#{$gtf_map_info_ex{$transcript_id}}){
       $count1 += $gtf_map_info_ex{$transcript_id}->[$index]->{"end"} - $gtf_map_info_ex{$transcript_id}->[$index]->{"start"} + 1;
     }
   }

   my $count2 = 0;
   open (G,"<$gtf_annotations") or die "Error[overlap_statistics.pl]! cannot read the $gtf_annotations file\n$!\n";
   foreach my $line (<G>){
     if ($line=~/^\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+/){
       $count2 += abs($2-$1) + 1;
     }
   }
   close G;

   system  "echo \"#hypereometric_test\" >> $outDir/${species}.ov.out";
   my $command_hyper = "perl $hypergeometric_test -all $allSize -count1 $count1 -count2 $count2 -jointcount $realOverlappingNts  >> $outDir/${species}.ov.out";
   (system "$command_hyper") == 0 or die "Error[overlap_statistics.pl]! cannot execute $command_hyper $!\n";
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





sub checkExonOverlap {
    my ($old , $new) = @_;
    my ($o_chr , $o_start , $o_end , $o_strand , $n_chr , $n_start , $n_end , $n_strand);
    if ($old=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+/){
	$o_chr    = $1;
	$o_start  = $2;
	$o_end    = $3;
	$o_strand = $4;
    }
    if ($new=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+/){
	$n_chr    = $1;
	$n_start  = $2;
	$n_end    = $3;
	$n_strand = $4;
    }
    my $overlap = 0;
    if (($o_chr eq $n_chr) && ($o_strand eq $n_strand)){
	$overlap = verifyOverlap ($o_start , $o_end , $n_start , $n_end);
    }

    return $overlap;
}
sub takeBlocks {
    my ($current_gtf) = @_;
    my $cmd_sort = 'sort -k 1,1 -k 7,7 -k 4,4n  ' . "$current_gtf";
    my @sorted_gtf = `$cmd_sort`; if ($?) {die "Error[reads_statistics.pl]! cannot run $cmd_sort\n";}
    my ($oldChr , $oldStart , $oldEnd , $oldStrand , %extractedTX , $foundExons , $exonsToBeFound , @current_block  );
    #initialize
    my $old_l = shift(@sorted_gtf);
    my $c = 0;
    if ($old_l =~/hitName \"([^\"]+)\"/){ 
	my $tx_id = $1;
	$extractedTX{$tx_id} = 1;
	push (@current_block , $tx_id);
	$exonsToBeFound += scalar(@{$gtf_map_info_ex{$tx_id}});
	$foundExons++;
    }
    my $overlap = 0;

    #LOOP
    while (my $l = shift(@sorted_gtf)){
	my $tx_id;
	if ($l =~/hitName \"([^\"]+)\"/){
	    $tx_id = $1;
	}

	if ((checkExonOverlap($l,$old_l) == 1) or ($foundExons < $exonsToBeFound)){
	    if (! defined $extractedTX{$tx_id}){ 
		$extractedTX{$tx_id} = 1;
		push (@current_block , $tx_id);
		$exonsToBeFound += scalar(@{$gtf_map_info_ex{$tx_id}});
	    }
	    $foundExons++;
	    $old_l      = $l;
	    next;
	}
	#save and reinitialize
	push (@{$blocks[$c]} , @current_block );
	$c++;
	@current_block             = ();
	%extractedTX               = ();
	$exonsToBeFound            = 0;
	push (@current_block , $tx_id);
	$exonsToBeFound += scalar(@{$gtf_map_info_ex{$tx_id}});
	$extractedTX{$tx_id} = 1;
	$foundExons          = 1;
	$old_l               = $l;
    }
    push (@{$blocks[$c]} , @current_block );
}
sub doBlocks {
    my $name = fileNameGenerator ('map.block.gtf');
    open (O,">$name") or die "Error[reads_statistics]! cannot create $name\n$!\n";
    foreach my $b (0..$#blocks){
	#print "block $b ";
	my (@allCoordinates , $strand , $chr , $outLine);
	foreach my $tx_id (@{$blocks[$b]}){
	    push (@allCoordinates , $gtf_map_info_tx{$tx_id}->{'start'});
	    push (@allCoordinates , $gtf_map_info_tx{$tx_id}->{'end'});
	    $strand = $gtf_map_info_tx{$tx_id}->{'strand'};
	    $chr    = $gtf_map_info_tx{$tx_id}->{'chr'};
	}
	my $blockStart  = min (@allCoordinates);
	my $blockEnd    = max (@allCoordinates);

	$outLine .= $chr              . "\t";
	$outLine .= 'read_statistics' . "\t";
	$outLine .= "block"           . "\t";
	$outLine .= $blockStart       . "\t";
	$outLine .= $blockEnd         . "\t";
	$outLine .= '.'               . "\t";
	$outLine .= $strand           . "\t";
	$outLine .= '.'               . "\t";
	$outLine .= "block \"$b\";";
	print O "$outLine\n";

	#take block internal deltas
	foreach my $tx_id (@{$blocks[$b]}){
	    foreach my $ex (@{$gtf_map_info_ex{$tx_id}}){
		my %hash;
		$hash{'delta_start'} = $ex->{'start'} -  $blockStart;
		$hash{'delta_end'}   = $ex->{'end'} -  $blockStart;
		push (@{$deltas4block[$b]{$tx_id}} , \%hash);
	    }
	}

	#save block start
	$blockStarts[$b] = $blockStart;
	$block_sizes[$b] = ($blockEnd - $blockStart) + 1 ;
    }
    close O;
    return $name;

}
sub blockGtf2exonGtf {
    my ($blockFile , $outFile) = @_;
    open (B,"<$blockFile") or die "Error[reads_statistics.pl]! cannot open $blockFile $!\n";
    open (O,">$outFile") or die "Error[reads_statistics.pl]! cannot open $outFile $!\n";

    while (my $current_block = <B>){
	#take shuffled block coordinates
	my ($p_chr , $p_start , $p_strand , $block_in);
	if ($current_block =~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+block \"([^\"]+)\";/){
	    $p_chr    = $1;
	    $p_start  = $2;
	    $p_strand = $3;
	    $block_in = $4;
	}

	#to do
	my $projection_delta = $p_start - $blockStarts[$block_in];

	#take internal transcript deltas
	foreach my $tx_id (@{$blocks[$block_in]}){
	    foreach my $ex_index (0..$#{$deltas4block[$block_in]{$tx_id}}){
		my $new_chr   = $p_chr;
		my $source    = $gtf_map_info_ex{$tx_id}[$ex_index]->{'source'};
		my $feature   = 'exon';
		my $new_start = $deltas4block[$block_in]{$tx_id}[$ex_index]{'delta_start'} + $p_start;
		my $new_end   = $deltas4block[$block_in]{$tx_id}[$ex_index]{'delta_end'}   + $p_start;
		my $score     = $gtf_map_info_ex{$tx_id}[$ex_index]->{'score'};
		my $strand    = $p_strand;
		my $frame     = $gtf_map_info_ex{$tx_id}[$ex_index]->{'frame'};
		my $group     = $gtf_map_info_ex{$tx_id}[$ex_index]->{'group'};
		print O "$new_chr\t$source\t$feature\t$new_start\t$new_end\t$score\t$strand\t$frame\t$group\n";
	    }

	}

    }
    close B;
    close O;
}
sub createRandomWithShuffleBed {
  my ($ran)       = @_;
  my $tmp = fileNameGenerator("_tmp_");
  my $cmd = "$shuffleBed -i $gtf_map_block -g $shuffleBed_chrSize";
  $cmd .= " -excl $shuffleBed_gap " if (defined $shuffleBed_gap);
  $cmd .= " > $tmp" ;
  (system "$cmd") == 0 or die "Error[overlap_statistics.pl]! cannot run $cmd \n$!\n";            # system "cat $tmp";
  blockGtf2exonGtf($tmp , $ran);
  system "rm $tmp";
}
sub createRandomGtf {
  my ($ran)   = @_;
  my $tmp = fileNameGenerator("_tmp_");
  my %mappedFeature;
  open (TMP,">>$tmp") or die "Error[overlap_statistics.pl]! cannot create $ran\n$!\n";
  foreach my $b (0..$#blocks){
    my $check  = 0;
    my ($r_chr , $r_start ,$r_end);
    while( $check  == 0 ){
      $r_chr   = random_chr();
      $r_start = random_start($r_chr , $block_sizes[$b]);
      $r_end   = $r_start +  $block_sizes[$b];
      if (%mappedFeature){
	$check   = ovCheck($r_chr , $r_start , $r_end , %mappedFeature);
      }
      else{$check = 1;}
    }
    my %hash = ('start' => $r_start , 'end' => $r_end);
    push (@{$mappedFeature{$r_chr}} , \%hash);
    print TMP "$r_chr\trandom\texon\t$r_start\t$r_end\t\.\t\.\t\.\tblock \"$b\"\n";
  }
  close TMP;
  blockGtf2exonGtf($tmp , $ran);
  system "rm $tmp";
}








##########EXTRAS

#OLD FUNCTION
# sub read_shuffleBed_chrSize {
#   open (RSC,"<$shuffleBed_chrSize") or die "Error! Cannot read $shuffleBed_chrSize";
#   foreach my $line (<RSC>){
#     if ($line=~/^(\S+)\s+(\S+)/){
#       $infoAssembly{$1} = $2;
#     }
#   }
# close RSC;
# }

#sub check_gtf_map {
#  my (%a , %b);
#  if ($mode eq 'exons'){
#    return $gtf_map ;
#  }
#  elsif ($mode eq 'transcripts'){
#    %a = readingGTFexonsLOCAL($gtf_map);
#    %b = gtfExons2Transcripts(%a);

#    $tx_gtf = fileNameGenerator("${outDir}/${species}.tx.gtf");
#    open (T,">$tx_gtf") or die "Error[overlap_statistics.pl]! Cannot create $tx_gtf\n$!\n";
#    foreach my $id (keys %b){
#      my $outLine;
#      $outLine .= $b{$id}->{'chr'}     . "\t";
#      $outLine .= $b{$id}->{'source'}  . "\t";
#      $outLine .= "transcript" . "\t";
#      $outLine .= $b{$id}->{'start'}   . "\t";
#      $outLine .= $b{$id}->{'end'}     . "\t";
#      $outLine .= $b{$id}->{'score'}   . "\t";
#      $outLine .= $b{$id}->{'strand'}  . "\t";
#      $outLine .= $b{$id}->{'frame'}   . "\t";
#      $outLine .= $b{$id}->{'group'}   . "\t";
#      print T "$outLine\n";
#    }
#    close T;
#  }
#  else {die "Error[overlap_statistics.pl]!  The mode $mode is not supported \n$!\n";}
#  undef (%a);
#  undef (%b);
#  return $tx_gtf;
#}
# sub reading_ex_gtf_annotations{
#   my (%gtf_annotations_intervalSizes);
#   open (G2,"<$ex_gtf_annotations") or die "Error[overlap_statistics.pl]! cannot read the $ex_gtf_annotations file\n";
#   my $feature = 0;
#   while (my $line = <G2>){
#     if ($line=~/^\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+/){
#       my $size = abs ($2 - $1);
#       $gtf_annotations_intervalSizes{$feature} = $size;
#       $feature++;
#     }
#   }
#   close G2;
#   return %gtf_annotations_intervalSizes;
# }


# sub doTranscriptGtf {
#     my (%gtf_map_info_tx) = @_;
#     my $name = fileNameGenerator ('map.tx.gtf');
#     open (O,">$name") or die "Error[reads_statistics]! cannot create $name\n$!\n";
#     foreach my $id (keys %gtf_map_info_tx){
# 	my $outLine;
# 	$outLine .= $gtf_map_info_tx{$id}->{'chr'}     . "\t";
# 	$outLine .= $gtf_map_info_tx{$id}->{'source'}  . "\t";
# 	$outLine .= "transcript" . "\t";
# 	$outLine .= $gtf_map_info_tx{$id}->{'start'}   . "\t";
# 	$outLine .= $gtf_map_info_tx{$id}->{'end'}     . "\t";
# 	$outLine .= $gtf_map_info_tx{$id}->{'score'}   . "\t";
# 	$outLine .= $gtf_map_info_tx{$id}->{'strand'}  . "\t";
# 	$outLine .= $gtf_map_info_tx{$id}->{'frame'}   . "\t";
# 	$outLine .= $gtf_map_info_tx{$id}->{'group'}   . "\t";
# 	print O "$outLine\n";
#     }
#     close O;
#     return $name;
# }
# sub randomSortTxId {
#     my (%hash) = @_;
#     my (%out , @sort);
#     foreach my $id (keys %hash){
# 	my $rNum = rand($numberTranscripts);
# 	$out{$rNum} = $id;
#     }
#     foreach my $i (sort { $b <=> $a } keys %out){
# 	push (@sort,$out{$i})
#     }
#     return @sort;
# }
# sub createRandomWithShuffleBed {
#   my ($ran)       = @_;
#     my $tmp_cov   = 0;
#     my $switch    = 0;

#     while ($switch == 0){
#       my $tmp = fileNameGenerator("_tmp_");
#       my @sortedTxId  = randomSortTxId(%gtf_map_info_tx);
#       my $cmd;
#       if ($gtfToRandomize eq 'gtf_map'){
# 	$cmd = "$shuffleBed -i $gtf_map_block -g $shuffleBed_chrSize";
#       }
#       elsif ($gtfToRandomize eq 'gtf_annotations'){
# 	$cmd = "$shuffleBed -i $ex_gtf_annotations -g $shuffleBed_chrSize";
#       }
#       $cmd .= " -excl $shuffleBed_gap " if (defined $shuffleBed_gap);
#       $cmd .= " > $tmp" ;
#       (system "$cmd") == 0 or die "Error[overlap_statistics.pl]! cannot run $cmd \n$!\n";            # system "cat $tmp";
#       my %ran_ex = transcriptGtf2exonGtf($tmp , \%gtf_map_info_ex);
#       (system "rm $tmp") == 0 or die "Error[overlap_statistics.pl]! cannot remove $tmp \n$!\n";
#       open (R,">>$ran")  or die "Error[overlap_statistics.pl]! cannot removeopen $ran \n$!\n";
#       while (($tmp_cov < $real_genome_coverage)&&(my $randomTx = shift (@sortedTxId))){
# 	foreach my $e (@{$ran_ex{$randomTx}}){
# 	  print R $e . "\n";
# 	}	                                                                                    #system "cat $ran";
# 	$tmp_cov = takeCoverage($ran);                                                              #print "$tmp_cov  $real_genome_coverage\n";
#       }
#       close R;
#       $switch = 1 unless ($tmp_cov < $real_genome_coverage);
#     }
# }






#old naming scheme related to overlap_statisticalRelevance_V4.pl
#gtf1 = gtf_map
#gtf2 = gtf_annotations





#LOOP
# . /etc/profile
# #!/bin/sh
# for VAR in  0 50 100 250 500 1000 1500 2000 5000
# do
# perl overlap_statisticalRelevance_V4.pl -gtf1 ~/Desktop/ncRNA-PROJECT/andrea/9000Lnc/gtf/gen3c.ex.processedTranscript_biggest.transcript.4each.gene.only_nonProteinOverlap.gtf -gtf2 ~/Desktop/ncRNA-PROJECT/vistaEnhancerParser/enhancers.gtf -assemblyDir ~/Desktop/allMammalInfo/human/assembly_hg19/ -gtfToRandomize gtf2 -extension_on_gtf2 $VAR  >> out
# done




#UCSC annotation example
#you can find an example here: ~/Desktop/ncRNA-PROJECT/vistaEnhancerParser/enhancer_annotations_4UCSC.gff
#bear in ming that the last field, the group, must be identical to the chromosome and that the "track name" line is compulsory

# browser position chr5:3186400-3187930
# track name=enhancerAnnotations description="vista_enhancer_annotations" color=0,0,255,
# chr5    vista   enhancer        3186439 3187926 0.00    .       .       chr5
# chr12   vista   enhancer        115124394       115125300       0.00    .       .       chr12
# chr4    vista   enhancer        105345575       105346895       0.00    .       .       chr4
# chr22   vista   enhancer        38394345        38395199        0.00    .       .       chr22
# chr2    vista   enhancer        162094895       162095451       0.00    .       .       chr2
# chr22   vista   enhancer        44398250        44399797        0.00    .       .       chr22
# chr15   vista   enhancer        58818261        58818883        0.00    .       .       chr15
# chr15   vista   enhancer        38159507        38161007        0.00    .       .       chr15
