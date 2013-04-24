#!/usr/bin/env perl
use strict;
use warnings;
use Cwd 'abs_path';
use Cwd 'cwd';
use File::Basename;
#Author: Giovanni Bussotti


#HELP
my $help  = 0;
my $pipelineDirName = cwd ();
my $scriptsDirName = dirname(abs_path($0));
my $allGenomeInfoDir = "${pipelineDirName}/allGenomeInfo";
unless (-d $allGenomeInfoDir) {(system "mkdir $allGenomeInfoDir") == 0 or die "Error[startPipeline.pl]! cannot create $allGenomeInfoDir  $!\n";}
my $experimentsDir = "${pipelineDirName}/experiments";
unless (-d $experimentsDir) {(system "mkdir $experimentsDir") == 0 or die "Error[startPipeline.pl]! cannot create $experimentsDir  $!\n";}
my $pipelineOutDir;
foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '-help') or ($ARGV[$field] eq '--help'));
}
if ($help > 0){
 help_message();
}
#TAKE OPTIONS
acceptedVariableSpace();
my ($repeatThreshold , $RNAalifold , $query_promoter_anchor , $splitExonerate , $ner , $splitGenome , $phylocsfName , $phyloCSFparametersName , $unstrandedGTFname , $querySplits4cluster4filtering , $uniprot , $link2referenceGenomeName , $overlapStrand , $cpus , $orf_score_threshold , $geneidParameterFile , $rfam , $annotationFile , $overlapDistance , $nr , $pfam , $blastx , $rpsblast , $codingPotential , $geneid , $referenceGenomeName , $query_gtf_name , $pre_processingName , $tcoffeeName , $xdformatName , $exonerateName , $genomesName , $chr_subseqName , $strategyName , $blastName , $queryName   , $exonerate_success_mode ,  $exonerate_lines_mode , $clusterName , $experimentName , $referencePhyloCsfName) = options();


#MAKE EMPTY DIRs
prepareDir();

#PRINT THE SCRIPT
open (P,">${pipelineOutDir}/RNAmapping_pipeline.txt") or die "Error[startPipeline.pl]! cannot create pipeline.log $!\n";
print P "#Pipeline started ";
my $date = `date`;
print P "$date";
print P "#commandline: $0 ";
foreach my $field (@ARGV){
  print P "$field ";
}
print P "\n";
$genomesName = abs_path($genomesName);
print P "NER=$ner\n";
print P "PHYLOCSFREFERENCE=$referencePhyloCsfName\n";
print P "SPLITGENOME=$splitGenome\n";
print P "PHYLOCSF=$phylocsfName\n";
print P "PHYLOCSFPARAMETERS=$phyloCSFparametersName\n";
print P "UNSTRANDEDGTF=$unstrandedGTFname\n";
print P "QUERYSPLITS4CLUSTER4FILTERING=$querySplits4cluster4filtering\n";
print P "UNIPROT=$uniprot\n";
print P "LINK2REFERENCEGENOMENAME=$link2referenceGenomeName\n";
print P "PREPROCESSING=$pre_processingName\n";
print P "GENOMES=$genomesName\n";
print P "QUERYFILE=$queryName\n";
print P "XDFORMAT=$xdformatName\n";
print P "EXONERATE=$exonerateName\n";
print P "CHR_SUBSEQ=$chr_subseqName\n";
print P "STRATEGY=$strategyName\n";
print P "BLAST=$blastName\n";
print P "PIPELINEDIR=$pipelineDirName\n";
print P "SCRIPTSDIR=${scriptsDirName}/scripts\n";
print P "EXONERATE_LINES_MODE=$exonerate_lines_mode\n";
print P "EXONERATE_SUCCESS_MODE=$exonerate_success_mode\n";
print P "CLUSTER=$clusterName\n";
print P "EXPERIMENT=$experimentName\n";
print P "TCOFFEE=$tcoffeeName\n";
print P "REFERENCEGENOME=$referenceGenomeName\n";
print P "QUERYGTF=$query_gtf_name\n";
print P "GENEID=$geneid\n";
print P "CODINGPOTENTIAL_CHECK=$codingPotential\n";
print P "RPSBLAST=$rpsblast\n";
print P "BLASTX=$blastx\n";
print P "PFAM=$pfam\n";
print P "NR=$nr\n";
print P "OVERLAPDISTANCE=$overlapDistance\n";
print P "OVERLAPSTRAND=$overlapStrand\n";
print P "ANNOTATION=$annotationFile\n";
print P "RFAM=$rfam\n";
print P "GENEIDPARAMETERFILE=$geneidParameterFile\n";
print P "ORFSCORETHRESHOLD=$orf_score_threshold\n";
print P "CPUS=$cpus\n";
print P "SPLITEXONERATE=$splitExonerate\n";
print P "QUERYPROMOTERANCHOR=$query_promoter_anchor\n";
print P "RNAALIFOLD=$RNAalifold\n";
print P "REPEATTHRESHOLD=$repeatThreshold\n";
print P "############\n";


print P "\n\n\n#PRE_PROCESSING\n";
print P 'if [ "$PREPROCESSING" = "on" ]; then $SCRIPTSDIR/createInput.pl -pipeline_dir $PIPELINEDIR -experiment $EXPERIMENT -reference_genome $REFERENCEGENOME -gtf $QUERYGTF -unstrandedGTF $UNSTRANDEDGTF -fasta $QUERYFILE -xdformat $XDFORMAT -blast $BLAST -blast_strategy $STRATEGY -cluster $CLUSTER  -splitGenome $SPLITGENOME -chr_subseq $CHR_SUBSEQ  -link2referenceGenome $LINK2REFERENCEGENOMENAME -exonerate_success_mode $EXONERATE_SUCCESS_MODE -exonerate_lines_mode $EXONERATE_LINES_MODE ; fi' . "\n";
print P 'CHECK=`echo $? |  tr -d "\n"`' . "\n";
print P 'if [ $CHECK -ne 0 ]; then exit ; fi' . "\n";
print P 'if [ "$PREPROCESSING" = "on" ]; then $SCRIPTSDIR/filtering.pl -geneid $GENEID -codingPotential_check $CODINGPOTENTIAL_CHECK -orf_score_threshold $ORFSCORETHRESHOLD -experiment $EXPERIMENT -pipeline_dir $PIPELINEDIR -rpsblast $RPSBLAST -blastx $BLASTX -pfam $PFAM -nr $NR -uniprot $UNIPROT -querySplits4cluster4filtering $QUERYSPLITS4CLUSTER4FILTERING -overlapDistance $OVERLAPDISTANCE -annotation $ANNOTATION -cluster $CLUSTER -rfam $RFAM -blast $BLAST -geneid_parameter $GENEIDPARAMETERFILE -overlapStrand $OVERLAPSTRAND ; fi' . "\n";
print P 'CHECK=`echo $? |  tr -d "\n"`' . "\n";
print P 'if [ $CHECK -ne 0 ]; then exit ; fi' . "\n";
print P 'if [ "$PREPROCESSING" = "on" ]; then echo "#The PRE_PROCESSING ran successfully. Now re-run the same command executePipeline.pl to do, if any, the further steps" ; exit ; fi' . "\n";



print P "\n\n\n#CREATE DATABASE\n";
print P '$SCRIPTSDIR/createDatabase.pl -genomes $GENOMES -xdformat $XDFORMAT -strategy $STRATEGY -pipeline_dir $PIPELINEDIR -splitGenome $SPLITGENOME -reference_genome $REFERENCEGENOME -experiment $EXPERIMENT'." \n";
print P 'CHECK=`echo $? |  tr -d "\n"`' . "\n";
print P 'if [ $CHECK -ne 0 ]; then exit ; fi' . "\n";



print P "\n\n\n#BLAST SEARCH\n";
print P 'rm $PIPELINEDIR/experiments/$EXPERIMENT/BLAST_OUT/* 2> /dev/null' . "\n";
print P '$SCRIPTSDIR/blastSearch.pl -query $QUERYFILE -blast $BLAST -blast_strategy $STRATEGY  -cluster $CLUSTER -splitGenome $SPLITGENOME -queryPromoterAnchor $QUERYPROMOTERANCHOR -experiment $EXPERIMENT -pipeline_dir $PIPELINEDIR' . "\n" ;
print P 'CHECK=`echo $? |  tr -d "\n"`' . "\n";
print P 'if [ $CHECK -ne 0 ]; then exit ; fi' . "\n";



print P "\n\n\n#EXONERATE REMAPPING\n";
print P 'echo "#exonerating..."' . "\n";
print P 'rm $PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/* 2> /dev/null' . "\n";
print P 'if [ $SPLITEXONERATE == "yes" ] && [ $CLUSTER == "on" ]; then bash $PIPELINEDIR/experiments/$EXPERIMENT/CONFIG/splitExonerate.sh ; for ID in `ls $PIPELINEDIR/allGenomeInfo`; do for SPLITMF2 in `ls $PIPELINEDIR/experiments/$EXPERIMENT/CLUSTER_FILES/splits_for_exonerate_$ID/grepOut_*`; do $SCRIPTSDIR/exonerateRemapping.pl -mf2 $SPLITMF2 -query $QUERYFILE -targetGenomeFolder $PIPELINEDIR/allGenomeInfo/$ID/chr/ -exonerate_lines_mode $EXONERATE_LINES_MODE -exonerate_success_mode $EXONERATE_SUCCESS_MODE -cluster $CLUSTER -experiment $EXPERIMENT -pipeline_dir $PIPELINEDIR -chr_subseq $CHR_SUBSEQ -ner $NER ;done ; done' . "\n";
print P 'CHECK=`echo $? |  tr -d "\n"`' . "\n";
print P 'else' . "\n";
print P 'for ID in `ls $PIPELINEDIR/allGenomeInfo`; do $SCRIPTSDIR/exonerateRemapping.pl -mf2 $PIPELINEDIR/experiments/$EXPERIMENT/BLAST_OUT/$ID.mf2 -query $QUERYFILE -targetGenomeFolder $PIPELINEDIR/allGenomeInfo/$ID/chr/ -exonerate_lines_mode $EXONERATE_LINES_MODE -exonerate_success_mode $EXONERATE_SUCCESS_MODE -cluster $CLUSTER -experiment $EXPERIMENT -pipeline_dir $PIPELINEDIR -chr_subseq $CHR_SUBSEQ -ner $NER ;done' . "\n";
print P 'CHECK=`echo $? |  tr -d "\n"`' . "\n";
print P 'fi' . "\n";
print P 'if [ $CHECK -ne 0 ]; then exit ; fi' . "\n";
print P "DONE_EXONERATE=`ls \$PIPELINEDIR/experiments/\$EXPERIMENT/CLUSTER_FILES/*___exonerate_done___  2> /dev/null | wc -l |  tr -d '\\n' ` \n";
print P 'if [ $SPLITEXONERATE == "yes" ] && [ $CLUSTER == "on" ]; then TO_BE_DONE=`ls $PIPELINEDIR/experiments/$EXPERIMENT/CLUSTER_FILES/splits_for_exonerate_*/x* | wc -l |  tr -d "\n" `' . "\n";
print P 'else' . "\n";
print P "TO_BE_DONE=`ls \$PIPELINEDIR/allGenomeInfo | wc -l |  tr -d '\\n' ` \n";
print P 'fi' . "\n";
print P 'while [ "$CLUSTER" = "on" -a "$DONE_EXONERATE" -lt "$TO_BE_DONE" ]; do sleep 20 ; echo "Waiting the exonerate job to exit the cluster. $DONE_EXONERATE already done..."; DONE_EXONERATE=`ls $PIPELINEDIR/experiments/$EXPERIMENT/CLUSTER_FILES/*___exonerate_done___  2> /dev/null | wc -l |  tr -d "\n" ` ; done; ';
print P 'if [ "$DONE_EXONERATE" -gt 0 ]; then rm  $PIPELINEDIR/experiments/$EXPERIMENT/CLUSTER_FILES/*___exonerate_done___ ; echo "#All Exonerate job exited the cluster"; fi' . "\n";
print P 'CHECK=`echo $? |  tr -d "\n"`' . "\n";
print P 'if [ $CHECK -ne 0 ]; then exit ; fi' . "\n";
print P 'if [ "$CLUSTER" = "on" ]; then for ERRORTEST in ERROR Error alloc Killed memory error; do CLUSTER_ERROR_CHECK=`grep $ERRORTEST $PIPELINEDIR/experiments/$EXPERIMENT/CLUSTER_FILES/*.ex.sh.e* | wc -l |  tr -d "\n"` ; if [ $CLUSTER_ERROR_CHECK -ne 0 ]; then echo "Error. A sanity check detected at least $CLUSTER_ERROR_CHECK times in $PIPELINEDIR/experiments/$EXPERIMENT/CLUSTER_FILES/ one of these possible error messages: ERROR Error error malloc Killed. Please detect it and rerun the script manually" ; exit ; fi ; done ; cat $PIPELINEDIR/experiments/$EXPERIMENT/CLUSTER_FILES/*.ex.sh.e* > $PIPELINEDIR/experiments/$EXPERIMENT/STDERR/exonerate; rm  $PIPELINEDIR/experiments/$EXPERIMENT/CLUSTER_FILES/*.ex.sh.e* ; fi' . "\n";
print P 'if [ $SPLITEXONERATE == "yes" ] && [ $CLUSTER == "on" ]; then for SPECIES in `ls $PIPELINEDIR/allGenomeInfo/`; do TOT=`ls $PIPELINEDIR/experiments/$EXPERIMENT/CLUSTER_FILES/splits_for_exonerate_$SPECIES/x* | wc -l | tr -d "\n"`; TOT=$((TOT-1)); X=-1; while [ $X -lt $TOT ] ; do X=$((X+1));  if [ $X -lt 10 ]; then NAMEex="$PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/grepOut_x0$X.$SPECIES.ex.gtf"; NAMEfa="$PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/grepOut_x0$X.$SPECIES.fa"; else NAMEex="$PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/grepOut_x$X.$SPECIES.ex.gtf"; NAMEfa="$PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/grepOut_x$X.$SPECIES.fa";fi;  cat $NAMEfa >> $PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/$SPECIES.fa ; cat $NAMEex >> $PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/$SPECIES.ex.gtf ; rm $NAMEfa ; rm $NAMEex ;done;done;fi' . "\n";
print P 'mkdir -p  $PIPELINEDIR/experiments/$EXPERIMENT/results' . "\n";
print P 'perl $SCRIPTSDIR/utility/repeatCoverage.pl $PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/ $PIPELINEDIR/experiments/$EXPERIMENT/results 20'  . "\n";
print P 'perl $SCRIPTSDIR/utility/repeatCoverage.pl $PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/ $PIPELINEDIR/experiments/$EXPERIMENT/results 80'  . "\n";
print P 'perl $SCRIPTSDIR/utility/repeatCoverage.pl $PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/ $PIPELINEDIR/experiments/$EXPERIMENT/results 100' . "\n";
print P 'perl $SCRIPTSDIR/utility/repeatCoverage.pl $PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/ $PIPELINEDIR/experiments/$EXPERIMENT/results $REPEATTHRESHOLD' . "\n";
print P 'for EX in `ls $PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/*.ex.gtf`; do NAME=`basename $EX | sed -e "s/.ex.gtf$/.all.rep100.ex.gtf/"| tr -d "\n"`; mv $EX $PIPELINEDIR/experiments/$EXPERIMENT/results/$NAME ; done' . "\n";
print P 'for FA in `ls $PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/*.fa`; do NAME=`basename $FA | sed -e "s/.fa$/.all.rep100.fa/"| tr -d "\n"`; mv $FA $PIPELINEDIR/experiments/$EXPERIMENT/results/$NAME ; done' . "\n";
print P 'for EX in `ls $PIPELINEDIR/experiments/$EXPERIMENT/results/*.ex.gtf`;do NAME=`echo $EX | sed -e "s/.ex.gtf$/.tx.gtf/" | tr -d "\n"`; perl $SCRIPTSDIR/utility/exonGTF_2_transcriptGTF.pl < $EX > $NAME;done ' . "\n";
print P 'for TX in `ls $PIPELINEDIR/experiments/$EXPERIMENT/results/*.tx.gtf`;do NAME=`echo $TX | sed -e "s/.tx.gtf$/.ge.gtf/" | tr -d "\n"`; perl $SCRIPTSDIR/utility/transcriptGTF_2_genesGTF.pl $TX > $NAME;done ' . "\n";
print P 'for EX in `ls $PIPELINEDIR/experiments/$EXPERIMENT/results/*.bh.rep$REPEATTHRESHOLD.ex.gtf`;do NAME=`basename $EX | sed -e "s/.bh.rep$REPEATTHRESHOLD.ex.gtf$/.ex.gtf/" | tr -d "\n"`; cp $EX $PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/$NAME ; done ' . "\n";
print P 'for FA in `ls $PIPELINEDIR/experiments/$EXPERIMENT/results/*.bh.rep$REPEATTHRESHOLD.fa`;do NAME=`basename $FA | sed -e "s/.bh.rep$REPEATTHRESHOLD.fa$/.fa/" | tr -d "\n"`; cp $FA $PIPELINEDIR/experiments/$EXPERIMENT/EXONERATE_OUT/$NAME ; done ' . "\n";




print P "\n\n\n#PREPARE_MFA\n";
print P 'rm $PIPELINEDIR/experiments/$EXPERIMENT/multifasta4EachTx/* 2> /dev/null' . "\n";
print P '$SCRIPTSDIR/prepare_mfa4alignment.pl -query $QUERYFILE  -experiment $EXPERIMENT -pipeline_dir $PIPELINEDIR ' . "\n";
print P 'CHECK=`echo $? |  tr -d "\n"`' . "\n";
print P 'if [ $CHECK -ne 0 ]; then exit ; fi' . "\n";



print P "\n\n\n#MULTIPLE SEQUENCE ALIGNMENT\n";
print P 'echo "#aligning..."' . "\n";
print P 'rm $PIPELINEDIR/experiments/$EXPERIMENT/outAlignments/* 2> /dev/null' . "\n";
print P 'ls $PIPELINEDIR/experiments/$EXPERIMENT/multifasta4EachTx/ | sed -e "s/.mfa$//" | xargs -I ID $TCOFFEE -in $PIPELINEDIR/experiments/$EXPERIMENT/multifasta4EachTx/ID.mfa -outfile $PIPELINEDIR/experiments/$EXPERIMENT/outAlignments/ID -method slow_pair -n_core $CPUS > /dev/null 2>> $PIPELINEDIR/experiments/$EXPERIMENT/STDERR/t_coffee ' ."\n";
print P 'rm -rf *.dnd' . "\n";
print P 'rm -rf $PIPELINEDIR/experiments/$EXPERIMENT/outAlignments/*.html' . "\n";

print P 'CHECK=`echo $? |  tr -d "\n"`' . "\n";
print P 'if [ $CHECK -ne 0 ]; then exit ; fi' . "\n";



print P "\n\n\n#EXTRACT SIMILARITY\n";
print P 'echo "#computing similarity..."' . "\n";
print P 'rm $PIPELINEDIR/experiments/$EXPERIMENT/outSim/* 2> /dev/null' . "\n";
print P 'ls $PIPELINEDIR/experiments/$EXPERIMENT/outAlignments/ | xargs -I ID sh -c "$TCOFFEE -other_pg seq_reformat -output sim -in $PIPELINEDIR/experiments/$EXPERIMENT/outAlignments/ID  > $PIPELINEDIR/experiments/$EXPERIMENT/outSim/ID" ' ."\n";
print P '$SCRIPTSDIR/sim2matrix.pl -query $QUERYFILE  -experiment $EXPERIMENT -pipeline_dir $PIPELINEDIR > $PIPELINEDIR/experiments/$EXPERIMENT/results/simMatrix.csv' . "\n";
print P 'R --slave --args $PIPELINEDIR/experiments/$EXPERIMENT/results/simMatrix.csv $PIPELINEDIR/experiments/$EXPERIMENT/results/ < $SCRIPTSDIR/utility/heatmap.R' . "\n";
print P 'mkdir -p  $PIPELINEDIR/experiments/$EXPERIMENT/consensusSecondaryStructures' . "\n";
print P 'for ALN_NAME in `ls $PIPELINEDIR/experiments/$EXPERIMENT/outAlignments/`; do $RNAALIFOLD -color -r -nc 0.5 -cv 0.6 -noLP $PIPELINEDIR/experiments/$EXPERIMENT/outAlignments/$ALN_NAME > $PIPELINEDIR/experiments/$EXPERIMENT/consensusSecondaryStructures/$ALN_NAME ; mv alirna.ps $PIPELINEDIR/experiments/$EXPERIMENT/consensusSecondaryStructures/$ALN_NAME.ps ;   done ' . "\n";
print P 'mkdir -p  $PIPELINEDIR/experiments/$EXPERIMENT/stockholmFormatAlignments' . "\n";
print P 'for ALN_NAME in `ls $PIPELINEDIR/experiments/$EXPERIMENT/outAlignments/`; do perl $SCRIPTSDIR/utility/viennaStructureAligner.pl $PIPELINEDIR/experiments/$EXPERIMENT/outAlignments/$ALN_NAME $PIPELINEDIR/experiments/$EXPERIMENT/consensusSecondaryStructures/$ALN_NAME > $PIPELINEDIR/experiments/$EXPERIMENT/stockholmFormatAlignments/$ALN_NAME ; done' . "\n";
print P 'mkdir -p  $PIPELINEDIR/experiments/$EXPERIMENT/compensatoryMutations' . "\n";
print P 'for STO_NAME in `ls $PIPELINEDIR/experiments/$EXPERIMENT/stockholmFormatAlignments/`; do $TCOFFEE -other_pg seq_reformat -in $PIPELINEDIR/experiments/$EXPERIMENT/stockholmFormatAlignments/$STO_NAME -action +alifold2analyze color_htm > $PIPELINEDIR/experiments/$EXPERIMENT/compensatoryMutations/$STO_NAME.html ; done' . "\n";
print P 'for STO_NAME in `ls $PIPELINEDIR/experiments/$EXPERIMENT/stockholmFormatAlignments/`; do echo "$STO_NAME " | tr -d "\n" >> $PIPELINEDIR/experiments/$EXPERIMENT/results/statistics_compensatory_mutations.txt ; $TCOFFEE -other_pg seq_reformat -in $PIPELINEDIR/experiments/$EXPERIMENT/stockholmFormatAlignments/$STO_NAME -action +alifold2analyze stat | grep @@ >> $PIPELINEDIR/experiments/$EXPERIMENT/results/statistics_compensatory_mutations.txt ; done' . "\n";


print P "\n\n\n#CSF\n";
print P 'echo "#computing the codon substitution frequency..."' . "\n";
print P 'rm $PIPELINEDIR/experiments/$EXPERIMENT/fasta_aln/* 2> /dev/null' . "\n";
print P '$SCRIPTSDIR/csf.pl -experiment $EXPERIMENT -pipeline_dir $PIPELINEDIR -t_coffee $TCOFFEE -phyloCSF $PHYLOCSF -phyloCSFparameters  $PHYLOCSFPARAMETERS -phyloCSFreference $PHYLOCSFREFERENCE' ."\n";



print P "\n\n\n" .'#END' . "\n";

close P;











###
#FUNCTIONS
sub options {
  my ($repeatThreshold , $RNAalifold , $query_promoter_anchor , $splitExonerate , $ner , $splitGenome , $phylocsfName , $phyloCSFparametersName , $unstrandedGTFname , $querySplits4cluster4filtering , $uniprot , $link2referenceGenomeName , $overlapStrand , $cpus , $orf_score_threshold , $geneidParameterFile , $rfam , $annotationFile , $overlapDistance , $nr , $pfam , $blastx , $rpsblast , $codingPotential , $geneid , $referenceGenomeName , $query_gtf_name , $pre_processingName , $tcoffeeName , $xdformatName , $exonerateName , $genomesName , $chr_subseqName , $strategyName , $blastName , $queryName   , $exonerate_success_mode ,  $exonerate_lines_mode  , $clusterName  , $experimentName , $referencePhyloCsfName);
  my $spyQuery_promoter_anchor         = 1;
  my $spyPhylocsfParameter             = 1;
  my $spyPhylocsf                      = 1;
  my $spyUnstrandedGTF                 = 1;
  my $spylink2referenceGenome          = 1;
  my $spyAnnotationFile                = 1;
  my $spyReferenceGenome               = 1;
  my $spyQuery_gtf_name                = 1;
  my $spyXdformat                      = 1;
  my $spyExonerate                     = 1;
  my $spyGenomes                       = 1;
  my $spyStrategy                      = 1;
  my $spyChr_subseq                    = 1;
  my $spyBlast                         = 1;
  my $spyQuery                         = 1;
  my $spyExperiment                    = 1;
  my $spyBlastconfig                   = 1;
  my $spyCluster                       = 1;
  my $spyClusterConfig                 = 1;
  my $spyTcoffee                       = 1;
  my $spyPre_processing                = 1;
  my $spyGeneid                        = 1;
  my $spyCodingPotential               = 1;
  my $spyRpsBlast                      = 1;
  my $spyBlastX                        = 1;
  my $spyPfam                          = 1;
  my $spyNr                            = 1;
  my $spyUniprot                       = 1;
  my $spyOverlapDistance               = 1;
  my $spyRfam                          = 1;
  my $spyGeneidParameterFile           = 1;
  my $spyExonerate_lines_mode          = 1;
  my $spyExonerate_success_mode        = 1;
  my $spyOrf_score_threshold           = 1;
  my $spyCpus                          = 1;
  my $spyOverlapStrand                 = 1;
  my $spyQuerySplits4cluster4filtering = 1;
  my $spyReferencePhyloCsf             = 1;
  my $spySplitGenome                   = 1;
  my $spyNer                           = 1;
  my $spySplitExonerate                = 1;
  my $spyRNAalifold                    = 1;
  my $spyRepeatThreshold               = 1;

  foreach my $field (0..$#ARGV){
    if ($ARGV[$field] eq '-repeatThreshold '){
	$repeatThreshold  = $ARGV[1+$field];
	$spyRepeatThreshold = 2;
	next;
    }
    if ($spyRepeatThreshold == 2){
	$spyRepeatThreshold = 3;
	next;
    }
    if ($ARGV[$field] eq '-RNAalifold'){
	$RNAalifold = $ARGV[1+$field];
	$spyRNAalifold = 2;
	next;
    }
    if ($spyRNAalifold == 2){
	$spyRNAalifold = 3;
	next;
    }
    if ($ARGV[$field] eq '-queryPromoterAnchor'){
	$query_promoter_anchor = $ARGV[1+$field];
	$spyQuery_promoter_anchor = 2;
	next;
    }
    if ($spyQuery_promoter_anchor == 2){
	$spyQuery_promoter_anchor = 3;
	next;
    }
    if ($ARGV[$field] eq '-splitExonerate'){
	$splitExonerate = $ARGV[1+$field];
	$spySplitExonerate = 2;
	next;
    }
    if ($spySplitExonerate == 2){
	$spySplitExonerate = 3;
	next;
    }
    if ($ARGV[$field] eq '-ner'){
	$ner = $ARGV[1+$field];
	$spyNer = 2;
	next;
    }
    if ($spyNer == 2){
	$spyNer = 3;
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
    if ($ARGV[$field] eq '-phyloCSFreference'){
	$referencePhyloCsfName = $ARGV[1+$field];
	$spyReferencePhyloCsf = 2;
	next;
    }
    if ($spyReferencePhyloCsf == 2){
	$spyReferencePhyloCsf = 3;
	next;
    }
    if ($ARGV[$field] eq '-phyloCSFparameters'){
	$phyloCSFparametersName = $ARGV[1+$field];
	$spyPhylocsfParameter = 2;
	next;
    }
    if ($spyPhylocsfParameter == 2){
	$spyPhylocsfParameter = 3;
	next;
    }
    if ($ARGV[$field] eq '-phyloCSF'){
	$phylocsfName = $ARGV[1+$field];
	$spyPhylocsf = 2;
	next;
    }
    if ($spyPhylocsf == 2){
	$spyPhylocsf = 3;
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
    if ($ARGV[$field] eq '-querySplits4cluster4filtering'){
	$querySplits4cluster4filtering = $ARGV[1+$field];
	$spyQuerySplits4cluster4filtering = 2;
	next;
    }
    if ($spyQuerySplits4cluster4filtering == 2){
      $spyQuerySplits4cluster4filtering = 3;
      next;
    }
   if ($ARGV[$field] eq '-uniprot'){
	$uniprot = $ARGV[1+$field];
	$spyUniprot = 2;
	next;
    }
    if ($spyUniprot == 2){
      $spyUniprot = 3;
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
    if ($ARGV[$field] eq '-overlapStrand'){
	$overlapStrand = $ARGV[1+$field];
	$spyOverlapStrand = 2;
	next;
    }
    if ($spyOverlapStrand == 2){
      $spyOverlapStrand = 3;
      next;
    }
    if ($ARGV[$field] eq '-cpus'){
	$cpus = $ARGV[1+$field];
	$spyCpus = 2;
	next;
    }
    if ($spyCpus == 2){
	$spyCpus = 3;
	next;
    }
    if ($ARGV[$field] eq '-orf_score_threshold'){
	$orf_score_threshold = $ARGV[1+$field];
	$spyOrf_score_threshold = 2;
	next;
    }
    if ($spyOrf_score_threshold == 2){
	$spyOrf_score_threshold = 3;
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
    if ($ARGV[$field] eq '-geneid_parameter'){
	$geneidParameterFile = $ARGV[1+$field];
	$spyGeneidParameterFile = 2;
	next;
    }
    if ($spyGeneidParameterFile == 2){
	$spyGeneidParameterFile = 3;
	next;
    }
    if ($ARGV[$field] eq '-rfam'){
      $rfam = $ARGV[1+$field];
      $spyRfam = 2;
      next;
    }
    if ($spyRfam == 2){
      $spyRfam = 3;
      next;
    }
    if ($ARGV[$field] eq '-annotation'){
	$annotationFile = $ARGV[1+$field];
	$spyAnnotationFile = 2;
	next;
    }
    if ($spyAnnotationFile == 2){
      $spyAnnotationFile = 3;
      next;
    }
    if ($ARGV[$field] eq '-overlapDistance'){
	$overlapDistance = $ARGV[1+$field];
	$spyOverlapDistance = 2;
	next;
    }
    if ($spyOverlapDistance == 2){
      $spyOverlapDistance = 3;
      next;
    }
    if ($ARGV[$field] eq '-nr'){
	$nr = $ARGV[1+$field];
	$spyNr = 2;
	next;
    }
    if ($spyNr == 2){
      $spyNr = 3;
      next;
    }
    if ($ARGV[$field] eq '-pfam'){
	$pfam = $ARGV[1+$field];
	$spyPfam = 2;
	next;
    }
    if ($spyPfam == 2){
      $spyPfam = 3;
      next;
    }
    if ($ARGV[$field] eq '-blastx'){
	$blastx = $ARGV[1+$field];
	$spyBlastX = 2;
	next;
    }
    if ($spyBlastX == 2){
      $spyBlastX = 3;
      next;
    }
    if ($ARGV[$field] eq '-rpsblast'){
	$rpsblast = $ARGV[1+$field];
	$spyRpsBlast = 2;
	next;
    }
    if ($spyRpsBlast == 2){
      $spyRpsBlast = 3;
      next;
    }
    if ($ARGV[$field] eq '-codingPotential_check'){
	$codingPotential = $ARGV[1+$field];
	$spyCodingPotential = 2;
	next;
    }
    if ($spyCodingPotential == 2){
	$spyCodingPotential = 3;
	next;
    }
    if ($ARGV[$field] eq '-geneid'){
	$geneid = $ARGV[1+$field];
	$spyGeneid = 2;
	next;
    }
    if ($spyGeneid == 2){
	$spyGeneid = 3;
	next;
    }
   if ($ARGV[$field] eq '-query_gtf'){
	$query_gtf_name = $ARGV[1+$field];
	$spyQuery_gtf_name = 2;
	next;
    }
    if ($spyQuery_gtf_name == 2){
	$spyQuery_gtf_name = 3;
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
   if ($ARGV[$field] eq '-pre_processing'){
	$pre_processingName = $ARGV[1+$field];
	$spyPre_processing = 2;
	next;
    }
    if ($spyPre_processing == 2){
	$spyPre_processing = 3;
	next;
    }
    if ($ARGV[$field] eq '-tcoffee'){
	$tcoffeeName = $ARGV[1+$field];
	$spyTcoffee = 2;
	next;
    }
    if ($spyTcoffee == 2){
	$spyTcoffee = 3;
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
    if ($ARGV[$field] eq '-cluster'){
	$clusterName = $ARGV[1+$field];
	$spyCluster = 2;
	next;
    }
    if ($spyCluster == 2){
	$spyCluster = 3;
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
    if ($ARGV[$field] eq '-blast'){
	$blastName = $ARGV[1+$field];
	$spyBlast = 2;
	next;
    }
    if ($spyBlast == 2){
	$spyBlast = 3;
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
    if ($ARGV[$field] eq 'exonerate'){
      $exonerateName = $ARGV[1+$field];
      $spyExonerate = 2;
      next;
    }
    if ($spyExonerate == 2){
      $spyExonerate = 3;
      next;
    }
    if ($ARGV[$field] eq 'chr_subseq'){
	$chr_subseqName = $ARGV[1+$field];
	$spyChr_subseq = 2;
	next;
    }
    if ($spyChr_subseq == 2){
      $spyChr_subseq = 3;
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
    if ($ARGV[$field] eq '-blast_strategy'){
	$strategyName = $ARGV[1+$field];
	$spyStrategy = 2;
	next;
    }
    if ($spyStrategy == 2){
	$spyStrategy = 3;
	next;
    }
  }


  if ((defined $overlapStrand) && ($overlapStrand ne 'yes') && ($overlapStrand ne 'no')){
    die "Error[startPipeline.pl]! The -overlapStrand parameter accepts either \"yes\" or \"no\"\n";
  }
  if ((defined $clusterName) && ($clusterName ne 'on') && ($clusterName ne 'off')){
    die "Error[startPipeline.pl]! The -cluster parameter accepts either \"on\" or \"off\"\n";
  }
  if ((defined $pre_processingName) && ($pre_processingName ne 'on') && ($pre_processingName ne 'off')){
    die "Error[startPipeline.pl]! The -pre_processing parameter accepts either \"on\" or \"off\"\n";
  }
  if ((defined $strategyName) && ($strategyName ne "wublastn") && ($strategyName ne "abblastn") && ($strategyName ne "wublastr") && ($strategyName ne "abblastr") && ($strategyName ne "wublastn_opt") && ($strategyName ne "abblastn_opt") ){
      die "Error[startPipeline.pl]! The supported blast flavors are wublastn abblastn wublastr abblastr wublastn_opt abblastn_opt\n";
  }
  if ((defined $codingPotential) && ($codingPotential ne 'on') && ($codingPotential ne 'off')){
    die "Error[startPipeline.pl]! The -codingPotential_check parameter accepts either \"on\" or \"off\"\n";
  }
  die "Error[startPipeline.pl]! You must provide the -genomes parameter. Please indicate either a folder containing the input genomes, either a file listing the pointers to the genomes files\n" if ($spyGenomes != 3);

  if (defined $queryName){
    $queryName   = abs_path($queryName);
    die "Error[startPipeline.pl]! The query file $queryName does not exist\n" unless (-f $queryName);
    my $check_header = `grep -c \">\" $queryName`;
    chomp $check_header;
    die "Error[startPipeline.pl]! The query file $queryName looks strange. It needs at least one header and one FASTA sequence\n" if ($check_header < 1);
    my $check_header_sort = `grep \">\" $queryName | sort | uniq | wc -l`;
    chomp $check_header_sort;
    die "Error[startPipeline.pl]! The query file $queryName contain duplicate IDs. Please amend it. The query identificators must be uniqe\n" if ($check_header != $check_header_sort);
  }
  if (defined $query_gtf_name){
    $queryName   = abs_path($queryName);
    die "Error[startPipeline.pl]! The query_gtf file $query_gtf_name does not exist\n" unless (-f $query_gtf_name);
    my $check_gtf = `grep -c  transcript_id $query_gtf_name`;
    chomp $check_gtf;
    die "Error[startPipeline.pl]! The query_gtf file $query_gtf_name file looks strange. It needs the transcript_id field\n" if ($check_gtf < 1);
  }
  if ((defined $query_promoter_anchor) && (defined $queryName)){
    my $hpa_tmp = fileNameGenerator("promoterAnchorHeaders");
    my $h_tmp   = fileNameGenerator("Headers");
    system ("grep \">\" $queryName | sort > $h_tmp");
    system ("grep \">\" $query_promoter_anchor | sort > $hpa_tmp");
    my $diffCheck = `diff $h_tmp $hpa_tmp`;
    if ($diffCheck){
      print "Error[startPipeline.pl]! the \"queries\" and the \"queries anchored by the promoter\" should have the same header!\n" ;
      system ("rm $h_tmp $hpa_tmp");
      exit;
    }
    system ("rm $h_tmp $hpa_tmp");
  }
  if (defined $genomesName){
    $genomesName = abs_path($genomesName);
    die "Error[startPipeline.pl]! The query file $genomesName does not exist\n" unless (-f $genomesName);
  }
  if ((defined $exonerate_lines_mode) && ($exonerate_lines_mode ne "exhaustive") && ($exonerate_lines_mode =~/\D/)){
      die "Error[exonerateRemapping.pl]! The supported exonerate_lines_mode modes are either \"exhaustive\" either an integer number of iteration\n";
  }
  if ((defined $exonerate_success_mode) && ($exonerate_success_mode ne "exhaustive") && ($exonerate_success_mode ne "ortholog") && ($exonerate_success_mode =~/\D/)){
      die "Error[exonerateRemapping.pl]! The supported exonerate_success_mode modes are either \"exhaustive\" either \"ortholog\" either  an integer number of iteration\n";
  }
  if (! defined $annotationFile){
    $overlapStrand = "no";
  }
  if ((! defined $overlapStrand) && (defined $annotationFile)){
    $overlapStrand = "yes";
  }
  if ((defined $phyloCSFparametersName) && ($phyloCSFparametersName ne '29mammals') && ($phyloCSFparametersName ne '12flies')){
    die "Error[startPipeline.pl]! You must provide the -phyloCSFparameters field with \'29mammals\' or \'12flies\'\n";
  }
  if ((defined $splitGenome) && ($splitGenome ne 'yes') && ($splitGenome ne 'no')){
    die "Error[startPipeline.pl]! -splitGenome parameter can be either \'yes\' or \'no\'\n";
  }
  if ((defined $ner)&&($ner ne 'no')&&($ner ne 'yes')){
    die "Error[startPipeline.pl]! -ner parameter can be either yes or no\n";
  }
  if ((defined $splitExonerate)&&($splitExonerate ne 'no')&&($splitExonerate ne 'yes')){
    die "Error[startPipeline.pl]! -splitExonerate parameter can be either yes or no\n";
  }

  $orf_score_threshold           = 20                                                    if (! defined $orf_score_threshold);
  $cpus                          = 1                                                     if (! defined $cpus);
  $exonerate_lines_mode          = 1000                                                  if (! defined $exonerate_lines_mode);
  $exonerate_success_mode        = 1                                                     if (! defined $exonerate_success_mode);
  $ner                           = 'no'                                                  if (! defined $ner);
  $annotationFile                = "none"                                                if (! defined $annotationFile);
  $overlapDistance               = 0                                                     if (! defined $overlapDistance);
  $nr                            = "none"                                                if (! defined $nr);
  $uniprot                       = "none"                                                if (! defined $uniprot);
  $querySplits4cluster4filtering = 10                                                    if (! defined $querySplits4cluster4filtering);
  $pfam                          = "none"                                                if (! defined $pfam);
  $blastx                        = 'blastx'                                              if (! defined $blastx);
  $rpsblast                      = 'rpsblast'                                            if (! defined $rpsblast);
  $codingPotential               = "off"                                                 if (! defined $codingPotential);
  $geneid                        = "geneid"                                              if (! defined $geneid);
  $queryName                     = 'none'                                                if (! defined $queryName);
  $query_gtf_name                = 'none'                                                if (! defined $query_gtf_name);
  $query_promoter_anchor         = 'none'                                                if (! defined $query_promoter_anchor);
  $splitGenome                   = 'no'                                                  if (! defined $splitGenome);
  $unstrandedGTFname             = 'off'                                                 if (! defined $unstrandedGTFname);
  $referenceGenomeName           = 'none'                                                if (! defined $referenceGenomeName);
  $pre_processingName            = 'off'                                                 if (! defined $pre_processingName);
  $tcoffeeName                   = 't_coffee'                                            if (! defined $tcoffeeName);
  $xdformatName                  = 'xdformat'                                            if (! defined $xdformatName);
  $exonerateName                 = 'exonerate'                                           if (! defined $exonerateName);
  $splitExonerate                = 'no'                                                  if (! defined $splitExonerate);
  $chr_subseqName                = 'chr_subseq'                                          if (! defined $chr_subseqName);
  $blastName                     = 'wu-blastn'                                           if (! defined $blastName);
  $strategyName                  = 'wublastn_opt'                                        if (! defined $strategyName);
  $rfam                          = 'none'                                                if (! defined $rfam);
  die "Error[startPipeline.pl]! Allowed cluster options are on|off \n"                   if ((defined $clusterName) && ($clusterName ne 'off') && ($clusterName ne 'on'));
  $clusterName                   = 'off'                                                 if (! defined $clusterName);
  $geneidParameterFile           = 'none'                                                if (! defined $geneidParameterFile);
  $link2referenceGenomeName      = 'off'                                                 if (! defined $link2referenceGenomeName);
  $phylocsfName                  = "PhyloCSF"                                            if (! defined $phylocsfName);
  $phyloCSFparametersName        = '29mammals'                                           if (! defined $phyloCSFparametersName);
  $referencePhyloCsfName         = "human"                                               if (! defined $referencePhyloCsfName);
  $RNAalifold                    = "RNAalifold"                                          if (! defined $RNAalifold);
  $repeatThreshold               = "20"                                                  if (! defined $repeatThreshold);



  die "Error[startPipeline.pl]! If you want to pre-process your query you need to provide the -reference_genome parameter\n" if (($pre_processingName eq 'on') && ($referenceGenomeName eq 'none'));
  sanity_chrHeaderSyntax($query_gtf_name , $referenceGenomeName) if (($query_gtf_name ne 'none') && ($referenceGenomeName ne 'none'));
  die "Error[startPipeline.pl]! You need to specify a query multi-FASTA file via -query \n" if (($pre_processingName eq 'off') && ($queryName eq 'none'));
  die "Error[startPipeline.pl]! You need to specify the -query and the -query_gtf parameters (or at least one of the two)\n" if (($pre_processingName eq 'on') && ($queryName eq 'none') && ($query_gtf_name eq 'none'));
  die "Error[startPipeline.pl]! If you want to use geneid you must specify with -geneid_parameter the parameter file estimated on your reference species. You can find a list of precomputed parameter files here: http://genome.crg.es/software/geneid/index.html#parameters either you can estimate a new one: http://genome.crg.es/software/geneid/training.html\n" if (($geneidParameterFile eq 'none') && ($codingPotential eq 'on'));
  die "Error[startPipeline.pl]! If you wanna enable the ortholog screening (by doing a reciprocal blast) you need to specify either the -query_gtf parameter either you have to run the pre_processing step (-pre_processing on). By doing this the script will be informed about the original query position on the reference genome\n" if (($exonerate_success_mode eq 'ortholog') && ($query_gtf_name eq 'none') && ($pre_processingName eq 'off'));
  die "Error[startPipeline.pl]! The -splitGenome function is meant to run on the cluster, please set -cluster on\n" if (($splitGenome eq 'yes') && ($clusterName ne 'on'));

  #DEFINE EXPERIMENT FOLDER
  if (! defined $experimentName){
    my $tmp_name_counter = 0;
    my $tmp_name;
    while (! $tmp_name || -d $tmp_name) {
      $tmp_name_counter++;
      $tmp_name = "${pipelineDirName}/experiments/exp_$tmp_name_counter";
    }
    $experimentName = "exp_$tmp_name_counter";
    (system "mkdir ${pipelineDirName}/experiments/exp_$tmp_name_counter") == 0 or die "Error[startPipeline.pl]! cannot create $experimentName $!\n" unless (-d "${pipelineDirName}/experiments/exp_$tmp_name_counter");
  }
  else{
    (system "mkdir ${pipelineDirName}/experiments/$experimentName") == 0 or die "Error[startPipeline.pl]! cannot create $experimentName $!\n" unless (-d "${pipelineDirName}/experiments/$experimentName");
  }
  $pipelineOutDir = "${pipelineDirName}/experiments/$experimentName";

  return ($repeatThreshold , $RNAalifold , $query_promoter_anchor , $splitExonerate , $ner , $splitGenome , $phylocsfName , $phyloCSFparametersName , $unstrandedGTFname , $querySplits4cluster4filtering , $uniprot , $link2referenceGenomeName , $overlapStrand , $cpus , $orf_score_threshold , $geneidParameterFile , $rfam , $annotationFile , $overlapDistance , $nr , $pfam , $blastx , $rpsblast , $codingPotential , $geneid , $referenceGenomeName , $query_gtf_name , $pre_processingName , $tcoffeeName , $xdformatName , $exonerateName , $genomesName , $chr_subseqName , $strategyName , $blastName , $queryName   , $exonerate_success_mode ,  $exonerate_lines_mode , $clusterName  , $experimentName , $referencePhyloCsfName);
}


sub acceptedVariableSpace {
  my %space = ('-repeatThreshold' => 1 , '-RNAalifold' => 1 , '-queryPromoterAnchor' => 1 , '-splitExonerate' => 1 , '-ner' => 1 , '-splitGenome' => 1 , '-phyloCSFparameters' => 1 , '-phyloCSF' => 1  , '-unstrandedGTF' => 1 , '-querySplits4cluster4filtering' => 1 ,'-uniprot' => 1 , '-link2referenceGenome' => 1 , '-overlapStrand' => 1  , '-cpus' => 1 , '-geneid_parameter' => 1 , '-rfam' => 1 , '-annotation' => 1 , '-overlapDistance' => 1 , '-nr' => 1 , '-pfam' => 1 , '-codingPotential_check' => 1 , '-blastx' => 1 , '-rpsblast' => 1 , '-query_gtf' => 1 , '-geneid' => 1  , '-reference_genome' => 1  , '-pre_processing' => 1  , '-tcoffee' => 1  , '-experiment' => 1  , '-cluster' => 1  ,'-exonerate_lines_mode' => 1  , '-exonerate_success_mode' => 1 , '-query' => 1  , '-chr_subseq' => 1   , '-blast' => 1   , '-xdformat' => 1   , '-exonerate' => 1   , '-blast_strategy' => 1   , '-genomes' => 1 , '-orf_score_threshold' => 1 , '-phyloCSFreference' => 1);
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[startPipeline.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}
sub sanity_chrHeaderSyntax {
  my ($query_gtf_name , $referenceGenomeName) = @_;
  print STDERR "#checking if query_gtf and reference genome have compatible chromosome names. This step could take time depenting on the reference genome size...\n";
  open (QGF,"<$query_gtf_name") or die "Error[startPipeline.pl]! cannot open $query_gtf_name $!\n";
  my $gtf_chr   = 0;
  my $gtf_NoChr = 0;
  while (my $l = <QGF>){
    next if ($l=~/^\s*$/);
    if ($l=~/^chr/){
      $gtf_chr = 1;
    }
    else {
      $gtf_NoChr = 1;
    }
  }
  close QGF;
  die "Error[startPipeline.pl]! in the query_gtf annotation file some chromosomes are preceded by a \'chr\' label, some are not. They should all have the same format (with or without \'chr\' label)\n" if ($gtf_chr == $gtf_NoChr);

  open (RGN,"<$referenceGenomeName") or die "Error[startPipeline.pl]! cannot open $referenceGenomeName $!\n";
  my $ref_chr   = 0;
  my $ref_NoChr = 0;
  while (my $l = <RGN>){
    if ($l=~/^>/){
      if ($l=~/^>chr/){
	$ref_chr = 1;
      }
      else  {
	$ref_NoChr = 1;
      }
    }
  }
  close RGN;
  die "Error[startPipeline.pl]! in the reference file some chromosomes are preceded by a \'chr\' label, some are not. They should all have the same format (with or without \'chr\' label)\n" if ($ref_chr == $ref_NoChr);

  die "Error[startPipeline.pl]! both reference and query_gtf files should have chromosomes with the same chromosome naming scheme (preceded or not by the \'chr\' label)\n" if ($gtf_chr != $ref_chr);
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
sub help_message {
my $helpMessage = "\nNAME
startPipeline.pl - Input the query RNAs and the target genomes to start a new RNAmapping experiment\n
SYNOPSIS
./startPipeline.pl -genomes [-query -query_gtf  -pre_processing -geneid_parameter -rfam -annotation -overlapDistance -nr -pfam -codingPotential_check -blastx -rpsblast -geneid -reference_genome -tcoffee -experiment -cluster -exonerate_mode -chr_subseq -blast -xdformat -exonerate -blast_strategy]\n
DESCRIPTION
   * startPipeline.pl takes as input a file listing:    input_genomes_file(absolute path)   name
   * By default (-pre_processing off) it takes a multi-FASTA file (-query) and screen the specified.
   * The user can pre_process the input queries (-pre_processing on) by entering either a fasta file (-query) either a gtf file (-query_gtf) together with the genome corresponding to the queries (-reference_genome)
   * By doing in this way the user can hunt new ncRNA sequences by applying different filters (none by default)

   OPTIONS
-------------------------------------
-general parameters:
   * -query                                                                          <Default: none>
   * -queryPromoterAnchor                                                            <Default: none>
   * -experiment                                                                     <Default: exp_XXX>
   * -cluster [on|off]                                                               <Default: off>
   * -cpus                                                                           <Default: 1>
-------------------------------------
-pre_processing parameters:
   * -pre_processing                                                                 <Default: off>
   * -query_gtf                                                                      <Default: none>
   * -reference_genome                                                               <Default: none>
   * -link2referenceGenome                                                           <Default: off>
   * -unstrandedGTF                                                                  <Default: off>
-filtering parameters:
   #annotation filter
     * -annotation                                                                   <Default: none>
     * -overlapDistance                                                              <Default: 0>
     * -overlapStrand                                                                <Default: yes>
   #coding potential filter
     * -codingPotential_check [on|off]                                               <Default: off>
     * -geneid                                                                       <Default: geneid>
     * -geneid_parameter                                                             <Default: none>
     * -orf_score_threshold                                                          <Default: 20>
   #rfam filter
     * -rfam                                                                         <Default: none>
   #blastX filter
     * -uniprot                                                                      <Default: none>
     * -nr                                                                           <Default: none>
     * -blastx                                                                       <Default: blastx>
     * -querySplits4cluster4filtering                                                <Default: 10>
   #rpsBlast filter
     * -pfam                                                                         <Default: none>
     * -rpsblast                                                                     <Default: rpsblast>
-------------------------------------
-database parameters:
   * -splitGenome                                                                    <Default: no>
   * -xdformat                                                                       <Default: xdformat>
-mapping parameters:
   * -blast_strategy [wublastn|abblastn|wublastr|abblastr|wublastn_opt|abblastn_opt] <Default: wublastn_opt>
   * -blast                                                                          <Default: wu-blastn>
   * -exonerate                                                                      <Default: exonerate>
   * -chr_subseq                                                                     <Default: chr_subseq>
   * -exonerate_lines_mode [integer|exhaustive]                                      <Default: 1000>
   * -exonerate_success_mode [integer|exhaustive|ortholog]                           <Default: 1>
   * -repeatThreshold                                                                <Default: 20>
   * -splitExonerate                                                                 <Default: no>
   * -ner                                                                            <Default: no>
   * -tcoffee                                                                        <Default: tcoffee>
   * -RNAalifold                                                                     <Default: RNAalifold>
-CSF parameters:
   * -phyloCSF                                                                       <Default: PhyloCSF>
   * -phyloCSFparameters                                                             <Default: 29mammals>
   * -phyloCSFreference                                                              <Default: human>

TROUBLESHOOTING
   * Read the MANUAL.txt file for a detailed description of the OPTIONS.
   * Both AB-Blast and WU-Blast packages are provided with their own xdformat. If needed it is therefore important that the user properly specifies, with -xdformat option, the proper one.
";
print "$helpMessage\n\n\n";
exit;
}




sub prepareDir {
  my $configDir = "${pipelineOutDir}/CONFIG" ;
  unless (-d $configDir) {(system "mkdir $configDir") == 0 or die "Error[startPipeline.pl]! cannot create $configDir  $!\n";}
  open (BC,">$configDir/blastConfig") or die "Error[startPipeline.pl]! cannot create the ${configDir}/blastConfig file $!\n";
  print BC "W=7 M=5 N=-4 Q=20 R=10 " if (($strategyName eq "abblastn_opt") or ($strategyName eq "wublastn_opt"));
  print BC "-e 0.00001 -cpus $cpus -filter=seg -lcfilter ";
  close BC;
  open (EEF,">$configDir/exonerateExtensionFile") or die "Error[startPipeline.pl]! Cannot create the ${configDir}/exonerateExtensionFile $!\n";
  print EEF "__arbitrary_extension__20000";
  close EEF;
  open (NBC,">$configDir/ncbiBlastXConfig") or die "Error[startPipeline.pl]! cannot create the ${configDir}/ncbiBlastxConfig file $!\n";
  print NBC "-evalue 0.0000000001 -num_alignments 1 -num_descriptions 1 -max_target_seqs 1 -num_threads $cpus ";
  close NBC;
  open (RBC,">$configDir/ncbiRPSblastConfig") or die "Error[startPipeline.pl]! cannot create the ${configDir}/ncbiRPSblastConfig file $!\n";
  print RBC "-evalue 0.0000000001 -num_alignments 1 -num_descriptions 1 -max_target_seqs 1 -num_threads $cpus ";
  close RBC;
  open (BR,">$configDir/blast4rfamFiltering") or die "Error[startPipeline.pl]! cannot create the ${configDir}/blast4rfamFiltering file $!\n";
  print BR "W=7 M=5 N=-4 Q=20 R=10 " if (($strategyName eq "abblastn_opt") or ($strategyName eq "wublastn_opt"));
  print BR "-e 0.00001 -cpus $cpus -filter=seg -lcfilter";
  close BR;

  open (GIC,">$configDir/geneidConfig") or die "Error[startPipeline.pl]! cannot create the ${configDir}/geneidConfig file $!\n";
  print GIC "-soW"; #-P $gene_id_param
  close GIC;
  open (ESF,">$configDir/splitExonerate.sh") or die "Error[startPipeline.pl]! Cannot create the ${configDir}/splitExonerate.sh $!\n";
  print ESF '. /etc/profile' . "\n";
  print ESF '#!/bin/sh'      . "\n";
  print ESF 'ALLGENOME=' . "$allGenomeInfoDir" . "\n";
  print ESF 'LINES=500'      . "\n";
  if ($queryName ne "none"){
    print ESF 'INPUT='."$queryName" . "\n";
  }
  else {
    print ESF 'INPUT='."${pipelineOutDir}/pre_processing/inputInfo/___query___.fa" . "\n";
  }
  print ESF 'for SPECIES in `ls $ALLGENOME`' . "\n";
  print ESF 'do' . "\n";
  print ESF 'mkdir -p '."${pipelineOutDir}".'/CLUSTER_FILES/splits_for_exonerate_$SPECIES' . "\n";
  print ESF "cd ${pipelineOutDir}".'/CLUSTER_FILES/splits_for_exonerate_$SPECIES/' . "\n";
  print ESF 'grep ">" $INPUT | tr -d ">" >' . " ${pipelineOutDir}".'/CLUSTER_FILES/splits_for_exonerate_$SPECIES/all.id'   . "\n";
  print ESF 'split -d  -l $LINES' .  " ${pipelineOutDir}".'/CLUSTER_FILES/splits_for_exonerate_$SPECIES/all.id'            . "\n";
  print ESF 'MF2='."${pipelineOutDir}".'/BLAST_OUT/$SPECIES.mf2' . "\n";
  print ESF 'for X in `ls x*`; do touch script_$X ; echo "#!/bin/sh" >> script_$X ;    echo  "for Y in \`cat $X\`;do grep -P \"^\$Y\\s+\" $MF2 2> /dev/null >> grepOut_$X.$SPECIES.mf2 ;done" >>  script_$X  ; echo "touch ___done_grep_$X.$SPECIES.___"  >> script_$X ;done' . "\n";
  print ESF 'for SCRIPT in `ls script*`; do qsub -V -cwd -o $PWD -e $PWD $SCRIPT -S $SHELL ; done' . "\n";
  print ESF 'cd ..' . "\n";
  print ESF 'done'  . "\n";

  print ESF 'DONE_GREP=`ls '."${pipelineOutDir}".'/CLUSTER_FILES/splits_for_exonerate_*/___done_grep_*  2> /dev/null | wc -l |  tr -d "\n"`' . "\n";
  print ESF 'TO_BE_DONE=`ls '."${pipelineOutDir}".'/CLUSTER_FILES/splits_for_exonerate_*/x*  2> /dev/null | wc -l |  tr -d "\n"`' . "\n";
  print ESF 'while [ "$DONE_GREP" -lt "$TO_BE_DONE" ]; do sleep 20 ; echo "Waiting the grep job to exit the cluster. $DONE_GREP already done..."; DONE_GREP=`ls '."${pipelineOutDir}".'/CLUSTER_FILES/splits_for_exonerate_*/___done_grep_*  2> /dev/null | wc -l |  tr -d "\n" ` ; done; ' . "\n";
  print ESF 'if [ "$DONE_GREP" -gt 0 ]; then rm  '."${pipelineOutDir}".'/CLUSTER_FILES/splits_for_exonerate_*/___done_grep_*  ; echo "#All grep job exited the cluster"; fi' . "\n";
  close ESF;

  my $clusterDir = "${pipelineOutDir}/CLUSTER_FILES";
  unless (-d $clusterDir) {(system "mkdir $clusterDir") == 0 or die "Error[startPipeline.pl]! cannot create $clusterDir  $!\n";}
  open (CC,">$configDir/clusterConfig") or die "Error[startPipeline.pl]! Cannot create the ${configDir}/clusterConfig $!\n";
  print CC "#COMMAND\n";
  print CC "qsub -V -o /dev/null -e ${pipelineOutDir}/CLUSTER_FILES -r y -cwd".' -S $SHELL  ##SCRIPT##';
  print CC "\n\n#HEADER\n";
  print CC '. /etc/profile' . "\n";
  print CC '#!/bin/sh' . "\n";
  close CC;

  my $stderrDir = "${pipelineOutDir}/STDERR";
  unless (-d $stderrDir) {(system "mkdir $stderrDir") == 0 or die "Error[startPipeline.pl]! cannot create $stderrDir  $!\n";}
  my $blastOut = "${pipelineOutDir}/BLAST_OUT";
  unless (-d $blastOut) {(system "mkdir $blastOut") == 0 or die "Error[startPipeline.pl]! cannot create $blastOut  $!\n";}
  my $exonerateOut = "${pipelineOutDir}/EXONERATE_OUT";
  unless (-d $exonerateOut) {(system "mkdir $exonerateOut") == 0 or die "Error[startPipeline.pl]! cannot create $exonerateOut  $!\n";}
  my $prepareOut = "${pipelineOutDir}/multifasta4EachTx";
  unless (-d $prepareOut) {(system "mkdir $prepareOut") == 0 or die "Error[startPipeline.pl]! cannot create $prepareOut  $!\n";}
  my $alignmentOut = "${pipelineOutDir}/outAlignments";
  unless (-d $alignmentOut) {(system "mkdir $alignmentOut") == 0 or die "Error[startPipeline.pl]! cannot create $alignmentOut  $!\n";}
  my $outSim = "${pipelineOutDir}/outSim";
  unless (-d $outSim) {(system "mkdir $outSim") == 0 or die "Error[startPipeline.pl]! cannot create $outSim  $!\n";}
  my $pre_processingOut = "${pipelineOutDir}/pre_processing";
  unless (-d $pre_processingOut) {(system "mkdir $pre_processingOut") == 0 or die "Error[startPipeline.pl]! cannot create $pre_processingOut  $!\n";}
  my $fasta_aln = "${pipelineOutDir}/fasta_aln";
  unless (-d $fasta_aln) {(system "mkdir $fasta_aln") == 0 or die "Error[startPipeline.pl]! cannot create $fasta_aln  $!\n";}
}

