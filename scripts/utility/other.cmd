
#EXEMPLE TO GENERATE SIZE DISTRIBUTIONS OF BOTH QUERY AND TARGET GENOMES, MEASURING THE COVERAGE FOR EACH TRANSCRIPT AND PLOT THE COVERAGES INTO AN HISTOGRAM
perl ~/my_script/GTF_manipulator_package/gtf2geneSizeDistribution.pl -gtfFile ../../../input/gencode.v7.allnoncoding.ex.all.gtf  -mode transcript | perl -ne 'chomp $_;if ($_=~/^(\S+)\.\d+ (\S+)/){print "$1 $2\n";}else{print"$_\n";}' > human.mature.tx.dist
perl ~/my_script/GTF_manipulator_package/gtf2geneSizeDistribution.pl -gtfFile mouse.ex.gtf  -mode transcript | perl -ne 'chomp $_;if ($_=~/^(\S+)\.\d+ (\S+)/){print "$1 $2\n";}else{print"$_\n";}' > mouse.mature.tx.dist
awk {'print $1'} mouse.mature.tx.dist | xargs -I X sh -c 'grep X human.mature.tx.dist' > human.mature.tx.dist.selection

to generate the coverages I used gnumeric,
then to generate the histogram plot I used R (on the dist file)




#TO CLEAN A BIT THE EXONERATE GTF OUTPUT, AND ASSIGNING THE CORRECT GENE_ID
# to take the original human transcript_id and gene_id and fix the problem of the divverent GENECODE version IDs
cat mouse.all.rep100.ex.gtf | perl -ne 'chomp $_; if($_=~/query_id \"([^"]+)\"/){$id=$1;} if($_=~/hitName \"([^"]+)\"/){$hname=$1;}  if($_=~/^(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+)/){$a=$1;} $gOut=`grep $id /users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/input/gencode.v10.long_noncoding_RNAs.tx.gtf `;chomp $gOut; if($gOut=~/transcript_id \"([^\"]+)\"/){$tx_id=$1;}  if($gOut=~/gene_id \"([^\"]+)\"/){$ge_id=$1;}    print "${a}gene_id \"$ge_id\"; transcript_id \"$tx_id\"; hitName \"$hname\";\n"' > mouse.clean.all.rep100.ex.gtf &

cat mouse.tx.gtf | perl -ne 'chomp $_; if($_=~/query_id \"([^"]+)\"/){$id=$1;} if($_=~/hitName \"([^"]+)\"/){$hname=$1;} if($_=~/^(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+)/){$a=$1;} $gOut=`grep $id /users/cn/gbussotti/Desktop/ncRNA-PROJECT/thomas2011_v7/input/gencode.v7.allnoncoding.tx.all.gtf`;chomp $gOut; if($gOut=~/transcript_id \"([^\"]+)\"/){$tx_id=$1;}  if($gOut=~/gene_id \"([^\"]+)\"/){$ge_id=$1;}    print "${a}gene_id \"$ge_id\"; transcript_id \"$tx_id\"; hitName \"$hname\";\n"'  > mouse4fantom.tx.gtf






#TO TAKE THE COVERAGE AND THE SIMILARITY DISTRIBUTIONS BETWEEN QUERY AND TARGET lncRNA
perl ~/my_script/GTF_manipulator_package/gtf2geneSizeDistribution.pl -gtfFile  /users/cn/gbussotti/Desktop/ncRNA-PROJECT/thomas2011_v7/input/gencode.v7.allnoncoding.ex.all.gtf  -mode transcript | perl -ne 'chomp $_;if ($_=~/^(\S+)\.\d+ (\S+)/){print "$1 $2\n";}else{print"$_\n";}' > human.mature.tx.dist

perl ~/my_script/GTF_manipulator_package/gtf2geneSizeDistribution.pl -gtfFile ~/Desktop/strategies4mappingLncRNA/abblastnOpt/experiments/exp_1/EXONERATE_OUT/cow.ex.gtf  -mode transcript | perl -ne 'chomp $_;if ($_=~/^(\S+)\.\d+ (\S+)/){print "$1 $2\n";}else{print"$_\n";}' > cow.mature.tx.dist
cat cow.mature.tx.dist | perl -ne 'chomp $_;if($_=~/^(\S+)\s+(\S+)\s*/){$id=$1;$num=$2;}else{die "ERROR\n"} $gOut=`grep $id human.mature.tx.dist`;chomp $gOut; if ($gOut=~/^\S+\s+(\S+)\s*/){$den=$1;}else{die"ERROR2!\n"} $cov=($num/$den)*100 ;  printf ("%.2f",$cov); print"\n";' | sort -n > cow.cov.dist
cat cow.mature.tx.dist | perl -ne 'chomp $_;if($_=~/^(\S+)\s+\S+\s*/){$id=$1}else{die "ERROR\n"} $gOut=`grep $id /users/cn/gbussotti/Desktop/strategies4mappingLncRNA/abblastnOpt/experiments/exp_1/simMatrix.csv`;chomp $gOut; if ($gOut=~/^[^,]+,([^,]+),[^,]+,[^,]+/){$sim=$1;}else{die"ERROR2!\n"} print"$sim\n";' | sort -n > cow.sim.dist

#then create the various histograms
cat plotHistogramForAdistribution.R | R --slave --args cow.cov.dist 10
cat plotHistogramForAdistribution.R | R --slave --args cow.sim.dist 10





#This oneliner tell for each fasta of a multiFASTA file how much is the REPEAT COVERAGE
#you can use it to print an histogram showing which fraction of your pipeR mapping is fully covered by repeats (or to show which fraction of input query was fully covered by repeats)
#one could consider to remove the target fully covered by repeats for instance
perl -e 'open(F,"<human.fa");while(<F>){chomp $_;if($_=~/>(\S+)/){if(defined $fraction){print "$fraction\n";} print"$1\t";   $size=0;$rep=0;next;}while($_=~/(\S)/g){$s=$1;next if ($s eq " "); $size++;if(lc($s)eq($s)){$rep++;} } $fraction=($rep/$size)*100;}  print "$fraction\n"; close F;'




# COMPARE THE NUMBER OF EXONS in the query and in the target for each hit:   
#  OUTPUT:  queryName hitName exonsInTheQuery exonsInTheTarget Difference
perl -e 'open(T,"</users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR/experiments/mouseENCODE_one2many_freyhult/EXONERATE_OUT/mouse.ex.gtf");@t=<T>;close T;  open(Q,"</users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/input/gencode.v10.long_noncoding_RNAs_foundInMouse.ex.gtf");@q=<Q>;close Q;          foreach $l (@t){if($l=~/hitName \"([^\"]+)\"/){$hit=$1;}if($hit=~/^(\S+)_hit\d+$/){$name=$1;}    $t_count=0;  foreach $ll(@t){if($ll=~/$hit\"/){$t_count++;}}   $q_count=0; foreach $lll(@q){if($lll=~/$name\"/){$q_count++;}}   $diff = $q_count - $t_count; print "$name $hit $q_count $t_count $diff\n";                                               }              ' | sort -k2 | uniq > exonNumberConservation 


#TO TAKE the number of query succesfully mapped in a certain species (e.g. mouse) 
for X in `cat ~/Desktop/ncRNA-PROJECT/mouseENCODE2012/input/gencode.v10.long_noncoding_RNAs.ex.gtf | ge_tr_hi_ID.pl tr | sort | uniq`; do echo "$X "| tr -d "\n" >> mouse_mapping_count ; grep -c $X ~/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012/experiments/mouseENCODE_one2many_freyhult/results/mouse.all.rep100.fa >> mouse_mapping_count ;done &

#take the mature transcript length distribution
perl -e 'open(F,"<../results/mouse.all.rep20.fa"); $tot=0; while ($l=<F>) {chomp $l;$l=~s/\s*$//; if($l=~/>/){if($tot>0){print"$tot\n";} $tot=0; next;} $tot+=length($l);               } print "$tot\n" ;close F;'

#cg content
perl -e 'open(F,"<../../../../input/gencode.v10.long_noncoding_RNAs.fa"); $tot=0; while ($l=<F>) {chomp $l;$l=~s/\s*$//; if($l=~/>/){if($tot>0){$res=($cg/$tot)*100;print"$res\n"; } $cg=0;$tot=0;  next;}     while($l=~/(\S)/g){$s=lc($1);$tot++; if(($s eq "c")or($s eq "g")){$cg++;}}} $res=($cg/$tot)*100;print"$res\n"; ' > cg_gencode &

#exon distribution
for X in `cat ../../../../input/gencode.v10.long_noncoding_RNAs.ex.gtf | ge_tr_hi_ID.pl tr | sort | uniq`; do grep -c $X ../../../../input/gencode.v10.long_noncoding_RNAs.ex.gtf; done > exonNumberDist_gencode


###########COMMANDS TO TAKE THE SENSITIVITY OF A MAPPING EXPERIMENT
#FETCH THE GENE_ID FOR THE EXONERATE OUTPUT AND MAKE TX.GTF
cat ../EXONERATE_OUT/cow.ex.gtf | perl -ne 'chomp $_; if($_=~/query_id \"([^"]+)\"/){$id=$1;}  if($_=~/^(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+)/){$a=$1;} $gOut=`grep $id /users/cn/gbussotti/Desktop/ncRNA-PROJECT/thomas2011_v7/input/gencode.v7.allnoncoding.tx.all.gtf`;chomp $gOut; if($gOut=~/transcript_id \"([^\"]+)\"/){$tx_id=$1;}  if($gOut=~/gene_id \"([^\"]+)\"/){$ge_id=$1;}    print "${a}gene_id \"$ge_id\"; transcript_id \"$tx_id\";\n"'  > cow.CORRECT_GENEID.ex.gtf &
perl ~/my_script/GTF_manipulator_package/exonGTF_2_transcriptGTF.pl < cow.CORRECT_GENEID.ex.gtf > cow.CORRECT_GENEID.tx.gtf

#COUNT THE NUMBER OF DIFFERENT GENES THAT WAS POSSIBLE TO MAP
echo "gene number " | tr -d "\n"
cat cow.CORRECT_GENEID.tx.gtf | ge_tr_hi_ID.pl ge | sort | uniq | wc -l

#SUM OF THE MATURE TRANSCRIPT LENGTH (nt)
echo "sum of mature transcript length " | tr -d "\n"
perl ~/my_script/GTF_manipulator_package/gtf2geneSizeDistribution.pl -gtfFile cow.CORRECT_GENEID.ex.gtf -mode transcript | cut -f2 -d" " | perl -ne 'chomp $_; $tot=0; $tot +=$_; while(<STDIN>){$tot +=$_;}print "$tot\n";'

#TAKE REAL COVERAGE
echo "real nucleotide coverage " | tr -d "\n"
perl ~/my_script/GTF_manipulator_package/exonGTF_2_ntCoverage.pl -exonGtf cow.CORRECT_GENEID.ex.gtf

#COUNT NOT OVERLAPPING LOCI
echo "strict unique loci " | tr -d "\n"
perl ~/my_script/script_funzionanti_vari/overlapping_locus.pl cow.CORRECT_GENEID.tx.gtf | wc -l

