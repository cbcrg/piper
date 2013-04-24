options(expressions=500000) 
library(gplots)  
library(ape)
library('fastcluster')
require(gee)
require(nlme)
require(lattice)


#CHOOSE THE INPUT FILE
test<-read.csv("simMatrix.csv",header=TRUE,row.names=1)
#test<-read.csv("fake.csv",header=TRUE,row.names=1)
testMat<-as.matrix(test)

MyTree <- read.tree("pruned_tree_of_life")
MyTreeUltrametric<-chronopl(MyTree, lambda=0.1)
tree <- as.hclust(MyTreeUltrametric)
dend <- as.dendrogram(tree)


#OUTPUT
bmp("heatMap_forecedTree.bmp", width=1000, height=1000)
x <- heatmap.2(testMat, Colv=dend ,col=colorpanel(40, "blue", "green", "yellow"),scale="none",key=TRUE,keysize=0.7,symkey=FALSE,density.info="none",trace="none",cexRow=0.5,cexCol=1.1 , dendrogram="column" )
dev.off()







#NOTES

#This script is meant to generate a conservation map forcing the tree of life. To do that you need to provide a newick formatted tree in a file named "pruned_tree_of_life".
#The file looks like this:

#pruned_tree_life_18_species
#(((((((((human:0.00670,chimp:0.00667):0.00225,orangutan:0.01832):0.01434,macaque:0.00785):0.02197,marmoset:0.06613):0.05759,(((mouse:0.08451,rat:0.09159):0.19777,guinea_pig:0.22563):0.01015,rabbit:0.11423):0.01531):0.02059,((cow:0.06180,pig:0.07900):0.02017,(horse:0.10940,(cat:0.09861,(panda:0.05123,dog:0.05123):0.05123):0.04984):0.00622):0.01167):0.02366,elephant:0.08224):0.23473,opossum:0.12569):0.07166,platypus:0.45659):0.10950;


#you might start from a bigger tree, like this:

#Full tree of life from Ensembl
#((((((((((((((((((((((human:0.0067,chimp:0.006667):0.00225,gorilla_gorilla:0.008825):0.00968,orangutan:0.018318):0.01434,(macaque:0.007853,?papio_hamadryas:0.007637):0.029618):0.021965,marmoset:0.066131):0.05759,tarsius_syrichta:0.137823):0.011062,(microcebus_murinus:0.092749,otolemur_garnettii:0.129725):0.035463):0.015494,tupaia_belangeri:0.186203):0.004937,(((((mouse:0.084509,rat:0.091589):0.197773,dipodomys_ordii:0.211609):0.022992,guinea_pig:0.225629):0.01015,spermophilus_tridecemlineatus:0.148468):0.025746,(rabbit:0.114227,ochotona_princeps:0.201069):0.101463):0.015313):0.020593,((((vicugna_pacos:0.107275,(tursiops_truncatus:0.064688,(cow:0.061796,?ovis_aries:0.061796):0.061796):0.025153):0.0201675,pig:0.079):0.0201675,((horse:0.109397,(cat:0.098612,(panda:0.051229,dog:0.051229):0.051229):0.049845):0.006219,(myotis_lucifugus:0.14254,pteropus_vampyrus:0.113399):0.033706):0.004508):0.011671,(erinaceus_europaeus:0.221785,sorex_araneus:0.269562):0.056393):0.021227):0.023664,(((elephant:0.082242,procavia_capensis:0.155358):0.02699,echinops_telfairi:0.245936):0.049697,(dasypus_novemcinctus:0.116664,choloepus_hoffmanni:0.096357):0.053145):0.006717):0.234728,(opossum:0.125686,macropus_eugenii:0.122008):0.2151):0.071664,platypus:0.456592):0.109504,((((gallus_gallus:0.041384,meleagris_gallopavo:0.041384):0.041384,anas_platyrhynchos:0.082768):0.082768,taeniopygia_guttata:0.171542):0.199223,anolis_carolinensis:0.489241):0.105143):0.172371,xenopus_tropicalis:0.855573):0.311354,(((tetraodon_nigroviridis:0.224159,takifugu_rubripes:0.203847):0.195181,(gasterosteus_aculeatus:0.316413,oryzias_latipes:0.48197):0.05915):0.32564,danio_rerio:0.730752):0.147949):0.526688,?petromyzon_marinus:0.526688),(ciona_savignyi:0.8,ciona_intestinalis:0.8)cionidae:0.6)chordata:0.2,(?apis_mellifera:0.9,(((?aedes_aegypti:0.25,?culex_quinquefasciatus:0.25):0.25,?anopheles_gambiae:0.5)culicinae:0.2,drosophila_melanogaster:0.8)diptera:0.1)endopterygota:0.7)coelomata:0.1,caenorhabditis_elegans:1.7)bilateria:0.3,saccharomyces_cerevisiae:1.9)fungi_metazoa_group:0.3);



#Then you need to prune it using t_coffee:
#t_coffee -other_pg seq_reformat -in species_tree.tree -in2 prune -action +tree_prune > pruned_tree_of_life

#prune file
#>cat
#>chimp
#>cow
#>dog
#>elephant
#>guinea_pig
#>horse
#>human
#>macaque
#>marmoset
#>mouse
#>opossum
#>orangutan
#>panda
#>pig
#>platypus
#>rabbit
#>rat





#WARNING:
# The columns of your heatmap should already be sorted by the tree of life since the heatmap is not going to sort the column according to the tree. It just put the tree on the top of the figure
# but the input matrix should already be sorted according to the tree of life.
# if this is not the situation you can sort it using awk, for example:
#awk -F "," {' print $1  "," $9  "," $3  "," $14  "," $10  "," $11  "," $12  "," $19  "," $7  "," $18  "," $4  "," $16  "," $8  "," $2  "," $15  "," $5  "," $6  "," $13  "," $17 '} fake.csv







