. /etc/profile
#!/bin/bash
#This script can be used to assign to a piper lncRNAout.ex.gtf output file and a protein exon file (relative to the target species) a biotype according with the GENCODE definition.
#see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431492/

#Example Commandline
# bash biotype_check.sh -l experiments/mouseENCODE_one2many_freyhult/results/mouse.all.rep20.ex.gtf -p ../input/full_and_proteinCoding_annotations/mm65.proteinCoding.ex.gtf -t  ~/my_script/GTF_manipulator_package/exonGTF_2_transcriptGTF.pl -i ~/my_script/GTF_manipulator_package/exonGTF_2_intronGTF.pl

while getopts ':l:p:t:i:h' OPTION ; do
case $OPTION in
     l)    LN=$OPTARG;;
     p)    PR=$OPTARG;;
     t)    ex2tx=$OPTARG;;
     i)    ex2in=$OPTARG;;

     h)    echo "USAGE of this program:"
           echo "   -l    lncRNA gtf exon file"
           echo "   -p    protein gtf exon file"
	   echo "   -t    exon to transcript perl script"
	   echo "   -i    exon to intron perl script"
         exit 0;;
     \?)    echo "Unknown argument \"-$OPTARG\"."
         echo $HELP ; exit 1;;
     :)    echo "Option \"-$OPTARG\" needs an argument."
         echo $HELP ; exit 1;;
     *)    echo "Sorry, I made a mistake when programming this.
"\"$OPTION\";;
esac
done


#exon to transcripts
perl $ex2tx -tag hitName < $LN > tmp_ln_tx
perl $ex2tx              < $PR > tmp_pr_tx
perl $ex2in -tag hitName < $LN > tmp_ln_in


#Antisense RNAs: Locus that has at least one transcript that intersects any exon of a protein-coding locus on the opposite strand, or published evidence of antisense regulation of a coding gene.
overlap tmp_ln_tx $PR -st -1 -o ov_biotype
grep -v "ov_feat2: 0" ov_biotype | perl -ne 'if($_=~/transcript_id \"([^\"]+)\";/){$id=$1; print "$id antisense_RNA\n";}  ' > out_biotype
rm ov_biotype


#LincRNA: Locus is intergenic noncoding RNA.
overlap tmp_ln_tx tmp_pr_tx -o ov_biotype
grep "ov_feat2: 0" ov_biotype | perl -ne 'if($_=~/transcript_id \"([^\"]+)\";/){$id=$1; print "$id lincRNA\n";}  ' >> out_biotype
rm ov_biotype

#Sense overlapping: Locus contains a coding gene within an intron on the same strand.
overlap tmp_ln_in tmp_pr_tx -st 1 -i 2 -o ov_biotype
grep -v "i2_feat2: 0" ov_biotype | perl -ne 'if($_=~/transcript_id \"([^\"]+)\";/){$id=$1; print "$id sense_overlapping\n";}  ' >> out_biotype
rm ov_biotype


#Sense intronic: Locus resides within intron of a coding gene but does not intersect any exons on the same strand
overlap tmp_ln_tx tmp_pr_tx -st 1 -o ov_biotype
grep -v "ov_feat2: 0" ov_biotype | perl -ne '$_=~s/ ov_feat2: \d+//; print $_;' > ___tmp___biotype
rm ov_biotype
overlap ___tmp___biotype $PR -st 1 -o ov_biotype
grep "ov_feat2: 0"  ov_biotype | perl -ne 'if($_=~/transcript_id \"([^\"]+)\";/){$id=$1; print "$id sense_intronic\n";}  ' >> out_biotype
rm ov_biotype ___tmp___biotype


#Pseudogenes: lncRNA exons overlap protein exons on the same strand (NOT A GENCODE TYPE)
overlap $LN $PR -st 1 -o ov_biotype
grep -v "ov_feat2: 0" ov_biotype | perl -ne 'if($_=~/hitName \"([^\"]+)\";/){$id=$1; print "$id pseudogenes\n";}  ' >> out_biotype
rm ov_biotype

#Antisense Intronic: Locus resides within intron of a coding gene on the other strand and  does not intersect any protein exons (NOT A GENCODE TYPE)
overlap tmp_ln_tx tmp_pr_tx -st -1 -o ov_biotype
grep -v "ov_feat2: 0" ov_biotype | perl -ne '$_=~s/ ov_feat2: \d+//; print $_;' > ___tmp___biotype
rm ov_biotype
overlap ___tmp___biotype $PR -st -1 -o ov_biotype
grep "ov_feat2: 0"  ov_biotype | perl -ne 'if($_=~/transcript_id \"([^\"]+)\";/){$id=$1; print "$id antisense_intronic\n";}  ' >> out_biotype
rm ov_biotype ___tmp___biotype

#Antisense overlapping: Locus contains a coding gene within an intron on the other strand. (NOT A GENCODE TYPE)
overlap tmp_ln_in tmp_pr_tx -st -1 -i 2 -o ov_biotype
grep -v "i2_feat2: 0" ov_biotype | perl -ne 'if($_=~/transcript_id \"([^\"]+)\";/){$id=$1; print "$id antisense_overlapping\n";}  ' >> out_biotype
rm ov_biotype


#concatenated: overlapping locus, but not overlapping exons in any strand (NOT A GENCODE TYPE)                                                                                                                      
overlap tmp_ln_tx tmp_pr_tx -o ov_biotype                                                                                                                                                      
grep -v "ov_feat2: 0" ov_biotype | perl -ne 'if($_=~/transcript_id \"([^\"]+)\";/){$id=$1;print "$id\n";}' > ___tmp___overlappingLoci                                                          
rm ov_biotype                                                                                                                                                                                  
overlap $LN $PR  -o ov_biotype                                                                                                                                                                 
grep -v "ov_feat2: 0"  ov_biotype | perl -ne 'if($_=~/hitName \"([^\"]+)\";/){$id=$1; print "$id,,,\n";}  ' > ___tmp___overlappingExons                                                        
for X in `cat ___tmp___overlappingLoci`; do 
GOUT=`grep $X,,, ___tmp___overlappingExons`; 
if [ -z "$GOUT" ]; then 
    GOUT=`grep "$X " out_biotype`; 
    if [ -z "$GOUT" ]; then 
	X=`echo $X |tr -d "\n"` ; 
	echo "$X concatenated";
    fi;
fi; 
done >> out_biotype
rm ___tmp___overlappingLoci ___tmp___overlappingExons

#UNKNOWN Sanity check
for X in `cat $LN | perl -ne 'if($_=~/hitName \"([^\"]+)\";/){$id=$1; print "$id\n";}' `; do
X=`echo $X | tr -d "\n"`;
GOUT=`grep "$X " out_biotype`;
if [ -z "$GOUT" ]; then
    #echo "$X processed_transcript" >> out_biotype;
    echo "ERROR! found unknown situation";
    exit 1;
fi;
done



###########

#REMOVE doubles
sort out_biotype | uniq > out_biotype2
mv out_biotype2 out_biotype

#Processed Transcripts: that is transcripts belonging to multiple biotypes
rm -rf out2 2> /dev/null
for X in `cat $LN | perl -ne 'if($_=~/hitName \"([^\"]+)\";/){$id=$1; print "$id\n";}' | sort | uniq `; do
X=`echo $X | tr -d "\n"`; 
COUNT=`grep -c "$X " out_biotype | tr -d "\n"`;  
if [ "$COUNT" -eq 1 ]; then
    grep "$X " out_biotype >> out2 ;
#| perl -ne 'if($_!~/\S+\s+1\s+/){print $_;}';
else 
    echo "$X processed_transcript" >> out2 ;
fi;
done

sort out2 | uniq > out22
mv out22 out2





perl -e 'open(F,"<'$LN'"); while ($l=<F>){chomp $l; if($l=~/hitName \"([^\"]+)\";/){$x=$1;} $gout=`grep "$x " out2 | cut -f2 -d " "`; chomp $gout; $l .= " transcript_type \"$gout\";\n"; print $l;  }  close F;'
rm tmp_ln_tx tmp_pr_tx  tmp_ln_in out2 out_biotype ov_biotype
