#This script is meant to compare using a wilkoxon test two pipeR heatMaps accepted as input
# it first converts each heatmap in an histogram (the two vompared distributions hust have the same number of features, bins to be compared)
# then each bin is normalized by the total number of features in each heatmap. This is necessary because you might be comparing heatmap having different sizes (i.e. generated out of different query amount)
# then the script run the wilkoxon test on the 2 distributions so generated
# it also returns a plot of the 2 density curves
#        cat compareHeatmaps.R | R --slave --args simMatrix.csv1 simMatrix.csv2

Args <- commandArgs();                #Args[4] will contain the first parameter
h1 = (Args[4])
h2 = (Args[5])

heatmap1 <-read.table(h1 , header = TRUE , sep = ",")
heatmap2 <-read.table(h2 , header = TRUE , sep = ",")






#TAKE DATA
allHeat1=c()
for (X in names(heatmap1)) {if (X == "species") next;    allHeat1 = c(allHeat1 , heatmap1[,X]) ; }
allHeat2=c()
for (X in names(heatmap2)) {if (X == "species") next;    allHeat2 = c(allHeat2 , heatmap2[,X]) ; }

#TAKE DISTRIBUTION AND NORMALIZE
h1=hist(allHeat1,plot=FALSE,breaks = 100)
h2=hist(allHeat2,plot=FALSE,breaks = 100)
normalizedDist1=c()
for (X in h1$counts) {val = (X/length(allHeat1))*100 ;  normalizedDist1=c(val, normalizedDist1) }
normalizedDist2=c()
for (X in h2$counts) {val = (X/length(allHeat2))*100 ;  normalizedDist2=c(val, normalizedDist2) }

#TEST
wilcox.test(normalizedDist1 , normalizedDist2)


#DENSITY PLOT
jpeg("heatmapComparison_density.jpeg", width=800, height=800 , res=120)
d_allHeat1 = density(allHeat1)
d_allHeat2 = density(allHeat2)
if (max(d_allHeat1$y) >= max(d_allHeat2$y)) {
	plot (d_allHeat1,col="blue",xlim=c(0,100),xaxs="i")  #not to extend over xlim limit
	lines(d_allHeat2,col="red")
}
if (max(d_allHeat1$y) < max(d_allHeat2$y)) {
	plot (d_allHeat2,col="red",xlim=c(0,100),xaxs="i")
	lines(d_allHeat1,col="blue")
}
legend("topright",c("heatMap1","heatMap2"),fill=c("blue","red"))
dev.off()



