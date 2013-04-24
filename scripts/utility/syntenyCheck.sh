. /etc/profile
#!/bin/sh
#INFO
#This bash script is meant to ease the syntenyCheck, and it is suitable for the RNAmapping pipeline output.
#To make a synteny comparison between a pair of genomes the user must specify the genome annotations in gtf format (include just the genes)
#You can easily retrieve by using biomart the fields: GENEID chrXX start end
#then reformat it to gtf by doing       cat biomartFile | perl -ne 'chomp $_;  if ($_=~/Ensembl/){next;};       if ($_=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*/){print "chr$2\tspecies\tgene\t$3\t$4\t0.00\t\.\t\.\tgene_id \"$1\"; transcript_id \"$1\"\n";}  ' > out.gtf
#The other file to provide is also achievable from Biomart. This is a list of all the orthologs. See checkSynteny.pl for further details
#Edit the VARIABLES to run this bash script

#COMMON VARIABLES
FINDCLOSESTFEATURE=/users/cn/gbussotti/my_script/script_funzionanti_vari/findClosestFeature.pl
CHECKSYNTENY=/users/cn/gbussotti/my_script/script_funzionanti_vari/checkSynteny.pl
EX2TX=/users/cn/gbussotti/my_script/GTF_manipulator_package/exonGTF_2_transcriptGTF.pl
OUT=syntenyWorkingFolder
#EXPERIMENT SPECIFIC VARIABLES
QUERY_GTF=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/roryCollaborators/input/decreasing.ex.gtf
TARGET_GTF=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/roryCollaborators/experiments/decreasing_normalHomolog/EXONERATE_OUT/outExonerate_mouse.gtf
QUERY_ALLGENE_ANNOTATIONS=/users/cn/gbussotti/Desktop/allMammalInfo/human/annotations_hg19/GRCh37.p5.ensemblGenes64.gtf
TARGET_ALLGENE_ANNOTATIONS=/users/cn/gbussotti/Desktop/allMammalInfo/mouse/annotations_mm9/musMusculus_NCBIM37.ensemblGenes.gtf
ORTHOLOGS_BIOMART=/home/gbussotti/DATASETS/human_mouse_orthologs_biomart

############################################
############################################
############################################

mkdir -p $OUT 
#make tx.gtf
$EX2TX < $QUERY_GTF > $OUT/query.tx.gtf
$EX2TX < $TARGET_GTF > $OUT/target.tx.gtf
#make intervals
cat $OUT/query.tx.gtf | perl -ne 'if ($_=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+/){print "$1 $2 $3\n";}' > $OUT/query.int
cat $OUT/target.tx.gtf | perl -ne 'if ($_=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+/){print "$1 $2 $3\n";}' > $OUT/target.int
#detect close genes
perl $FINDCLOSESTFEATURE -allowOverlap -neighbors 3 -strategy neighbors -mode iVSg -queries $OUT/query.int -gtf $QUERY_ALLGENE_ANNOTATIONS > $OUT/query.close.genes
perl $FINDCLOSESTFEATURE -allowOverlap -neighbors 3 -strategy neighbors -mode iVSg -queries $OUT/target.int -gtf $TARGET_ALLGENE_ANNOTATIONS > $OUT/target.close.genes
#prepare idCorrespondenceFile
cat $OUT/target.tx.gtf | perl -ne 'if($_=~/transcript_id \"([^\"]+)\"/){$txID=$1;}if($_=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+/){$targetCODE="$1_$2_$3";}  $gOut=`grep $txID '$OUT/query.tx.gtf'`; if($gOut=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+/){$queryCODE="$1_$2_$3";}   print "$queryCODE $targetCODE\n";' > $OUT/idCorrespondenceFile
#prepare intervalCode TO txID
cat $OUT/target.tx.gtf | perl -ne 'if($_=~/transcript_id \"([^\"]+)\"/){$txID=$1;}if($_=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+/){$targetCODE="$1_$2_$3"; print "$txID $targetCODE\n";}' > $OUT/target.codes
cat $OUT/query.tx.gtf | perl -ne 'if($_=~/transcript_id \"([^\"]+)\"/){$txID=$1;}if($_=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+/){$queryCODE="$1_$2_$3"; print "$txID $queryCODE\n";}' > $OUT/query.codes 
#check synteny
echo "transcriptID queryInterval targetInterval numberOfOrthologs orthologPairsFound";
perl $CHECKSYNTENY -strategy id_match -mode gVSg -idCorrespondence $OUT/idCorrespondenceFile -hitGtf $TARGET_GTF -targetCloseGenes $OUT/target.close.genes -closeGenes $OUT/query.close.genes -allOrthologList $ORTHOLOGS_BIOMART | perl -ne 'if ($_=~/^(\S+)/){$queryCODE=$1;$gOut=`grep $queryCODE '$OUT/query.codes'`; if ($gOut=~/^(\S+)/){$txID=$1;print "$txID $_";}} ' 





