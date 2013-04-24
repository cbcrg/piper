. /etc/profile
#!/bin/bash
#This script can be used to compare for each transcript the query biotype with the target biotype
#Biotype info can be found here:
#http://www.gencodegenes.org/gencode_biotypes.html
#Both the query and the target files MUST have the transcript_type field
#As biotype can change with time, you need to provide also biotype equivalence file with this syntax:
#queryBiotype(TAB)targetBiotype,,,

# for instance:
#ambiguous_orf   pseudogenes,,,
#lincRNA         lincRNA,,,
#processed_transcript    processed_transcript,,,
#sense_intronic  sense_intronic,,,
#sense_overlapping       sense_overlapping,,,
#antisense               antisense_RNA,,,

#Example Commandline
#bash biotype_compare.sh -q ../input/gencode.v10.long_noncoding_RNAs.ex.gtf -t out_biotype_check.ex.gtf -e equivalences


while getopts ':q:t:e:h' OPTION ; do
case $OPTION in
     q)    Q=$OPTARG;;
     t)    T=$OPTARG;;
     e)    E=$OPTARG;;
     h)    echo "USAGE of this program:"
           echo "   -q    lncRNA gtf exon query file"
           echo "   -t    lncRNA gtf exon target file as returned by pipeR and converted by biotype_check.sh"
	   echo "   -e    transcript_type equivalences"
         exit 0;;
     \?)    echo "Unknown argument \"-$OPTARG\"."
         echo $HELP ; exit 1;;
     :)    echo "Option \"-$OPTARG\" needs an argument."
         echo $HELP ; exit 1;;
     *)    echo "Sorry, I made a mistake when programming this.
"\"$OPTION\";;
esac
done



rm -rf missingEquivalences  2> /dev/null
for X in `cat $T | perl -ne 'if($_=~/hitName \"([^\"]+)\"/){print "$1\n";}' | sort | uniq `; do 
    TT=`grep -m1 "$X\"" $T  | perl -ne 'if($_=~/transcript_type \"([^\"]+)\"/){print "$1";}'`;
    ID=`grep -m1 "$X\"" $T  | perl -ne 'if($_=~/transcript_id \"([^\"]+)\"/){print "$1";}'`;
    TTE=`grep $TT,,, $E | cut -f1`;

    QT=`grep  -m1 "$ID\"" $Q | perl -ne 'if($_=~/transcript_type \"([^\"]+)\"/){print "$1";}'`;
    CHECK=`grep -P "$QT\s+" $E`;
    if [ -z "$CHECK" ]; then echo $ID >> missingEquivalences; continue; fi;
    if [ "$TTE" == "$QT" ]; then echo "the target $X has the same biotype as the query $ID, $TTE"; else echo "the target $X has different biotype than the query $ID, $TTE"; fi;
done