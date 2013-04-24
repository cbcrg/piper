. /etc/profile
#!/bin/bash
#This simple bash script can be used to check the number of bh homologs found by pipeR in a certain species, check the predictions that are expressed, check if the promoter anchored mapping of the same method can rescue the mapping, and check if the promoter anchored predictions are expressed.

ALL=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/input/gencode.v10.long_noncoding_RNAs.fa
EXP_THRESHOLD=0.1

##BLASTN
DB1=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012/experiments/mouseENCODE_one2many_blastn/results/mouse.bh.rep20.ex.gtf
DB2=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012/experiments/mouseENCODE_one2many_blastn_anchor/results/mouse.bh.rep20.ex.gtf
EXP1=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012/experiments/mouseENCODE_one2many_blastn/qualityCheck/reads_statistics/bamCSHLmouse_VS_mouse_all_rep20/transcripts.rpkm
EXP2=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012/experiments/mouseENCODE_one2many_blastn_anchor/qualityCheck/reads_statistics/bamCSHLmouse_VS_mouse_all_rep20/transcripts.rpkm



echo "found by blastN";
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
    if [[ "$gout1" ]];
    then echo "$ID" ;
    fi;
done  | wc -l

echo "found by blastN and expressed over the rpkm thresold";
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
    if [[ "$gout1" ]];
    then EXPRESSION=`grep ${ID}_hit $EXP1  | sort -rnk2 | head -1 | cut -f2 | tr -d "\n"`  ;
	if [[ "$EXPRESSION" ]]; then 
	    perl -e 'if ('$EXPRESSION' >= '$EXP_THRESHOLD'){print "the highest expression for '$ID' is '$EXPRESSION'\n"}'; 
	    #CHECK=`echo $? |  tr -d "\n"`; 
            # if [ $CHECK -ne 0 ]; then echo $ID; echo $EXPRESSION ; exit ; fi    
	fi;         
    fi;
done  | wc -l  


echo "recovered in blastN by the promoter";
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
    if [ -z "$gout1" ]; 
    then gout2=`grep -m 1 ${ID}_ $DB2  | tr -d "\n"`;
	if [[ "$gout2" ]]; then echo "$ID" ; fi;
    fi;
done | wc -l

echo "recovered in blastN by the promoter and expressed over the rpkm thresold";
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
    if [ -z "$gout1" ]; 
    then gout2=`grep -m 1 ${ID}_hit $DB2  | tr -d "\n"`;
	if [[ "$gout2" ]]; 
	then EXPRESSION=`grep ${ID}_hit $EXP2  | sort -rnk2 | head -1 | cut -f2 | tr -d "\n"`  ;
	    if [[ "$EXPRESSION" ]]; then 
		perl -e 'if ('$EXPRESSION' >= '$EXP_THRESHOLD'){print "the highest expression for '$ID' if '$EXPRESSION'\n"}';
	    fi;
	fi;
    fi;
done | wc -l








##BLASTNOPT
DB1=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012/experiments/mouseENCODE_one2many_freyhult/results/mouse.bh.rep20.ex.gtf
DB2=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012/experiments/mouseENCODE_one2many_freyhult_anchor/results/mouse.bh.rep20.ex.gtf
EXP1=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012/experiments/mouseENCODE_one2many_freyhult/qualityCheck/reads_statistics/bamCSHLmouse_VS_mouse_all_rep20/transcripts.rpkm
EXP2=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012/experiments/mouseENCODE_one2many_freyhult_anchor/qualityCheck/reads_statistics/bamCSHLmouse_VS_mouse_all_rep20/transcripts.rpkm

echo "found by blastN opt";
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
if [[ "$gout1" ]];
    then echo "$ID" ;
fi;
done  | wc -l

echo "found by blastN opt and expressed over the rpkm thresold"
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
if [[ "$gout1" ]];
    then EXPRESSION=`grep ${ID}_hit $EXP1  | sort -rnk2 | head -1 | cut -f2 | tr -d "\n"`  ;
    if [[ "$EXPRESSION" ]]; then 
	perl -e 'if ('$EXPRESSION' >= '$EXP_THRESHOLD'){print "the highest expression for '$ID' if '$EXPRESSION'\n"}';
    fi;
fi;
done   | wc -l

echo "recovered in blastN opt by the promoter"
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
if [ -z "$gout1" ]; 
    then gout2=`grep -m 1 ${ID}_hit $DB2  | tr -d "\n"`;
    if [[ "$gout2" ]]; then echo "$ID" ; fi;
fi;
done | wc -l

echo "recovered in blastN opt by the promoter and expressed over the rpkm thresold";
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
    if [ -z "$gout1" ]; 
    then gout2=`grep -m 1 ${ID}_hit $DB2  | tr -d "\n"`;
	if [[ "$gout2" ]]; 
	then EXPRESSION=`grep ${ID}_hit $EXP2  | sort -rnk2 | head -1 | cut -f2 | tr -d "\n"`  ;
	    if [[ "$EXPRESSION" ]]; then 
		perl -e 'if ('$EXPRESSION' >= '$EXP_THRESHOLD'){print "the highest expression for '$ID' if '$EXPRESSION'\n"}';
	    fi;
	fi;
    fi;
done | wc -l


##BLASTR
DB1=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012_blastR/experiments/mouseENCODE_one2many_blastr/results/mouse.bh.rep20.ex.gtf
DB2=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012_blastR/experiments/mouseENCODE_one2many_blastrAnchored/results/mouse.bh.rep20.ex.gtf
EXP1=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012_blastR/experiments/mouseENCODE_one2many_blastr/qualityCheck/reads_statistics/bamCSHLmouse_VS_mouse_all_rep20/transcripts.rpkm
EXP2=/users/cn/gbussotti/Desktop/ncRNA-PROJECT/mouseENCODE2012/pipeR_3102012_blastR/experiments/mouseENCODE_one2many_blastrAnchored/qualityCheck/reads_statistics/bamCSHLmouse_VS_mouse_all_rep20/transcripts.rpkm

echo "found by blastR";
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
if [[ "$gout1" ]];
    then echo "$ID" ;
fi;
done  | wc -l

echo "found by blastR and expressed over the rpkm thresold"
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
if [[ "$gout1" ]];
    then EXPRESSION=`grep ${ID}_hit $EXP1  | sort -rnk2 | head -1 | cut -f2 | tr -d "\n"`  ;
    if [[ "$EXPRESSION" ]]; then 
	perl -e 'if ('$EXPRESSION' >= '$EXP_THRESHOLD'){print "the highest expression for '$ID' if '$EXPRESSION'\n"}';
    fi;
fi;
done   | wc -l

echo "recovered in blastR by the promoter"
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
if [ -z "$gout1" ]; 
    then gout2=`grep -m 1 ${ID}_hit $DB2  | tr -d "\n"`;
    if [[ "$gout2" ]]; then echo "$ID" ; fi;
fi;
done | wc -l

echo "recovered in blastR by the promoter and expressed over the rpkm thresold";
for ID in `grep ">" $ALL |sed -e "s/>//"`; do   gout1=`grep -m 1 ${ID}_hit $DB1 | tr -d "\n"`;
    if [ -z "$gout1" ]; 
    then gout2=`grep -m 1 ${ID}_hit $DB2  | tr -d "\n"`;
	if [[ "$gout2" ]]; 
	then EXPRESSION=`grep ${ID}_hit $EXP2  | sort -rnk2 | head -1 | cut -f2 | tr -d "\n"`  ;
	    if [[ "$EXPRESSION" ]]; then 
	    perl -e 'if ('$EXPRESSION' >= '$EXP_THRESHOLD'){print "the highest expression for '$ID' if '$EXPRESSION'\n"}';
	    fi;
	fi;
    fi;
done | wc -l



