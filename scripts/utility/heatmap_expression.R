options(expressions=500000) #this option is needed if you have many features in your matrix and you have to tell R that it is not an infinite loop, but just a bit longer
###example commandLine: R --slave --args simMatrix.csv experiments/exp_1/results/ < qualityCheck/utility/heatmap.R
Args <- commandArgs();
simMatrix.csv=Args[4]
outDir=Args[5]


#LIBRARY NEEDED IF THE HEATMAP IS HUGE (i.e. 17000 transcripts)
# if the package is not installed use the command  install.packages("fastcluster") to install it
library('fastcluster')



#CHOOSE THE INPUT FILE
#test<-read.csv("heatMap_noAllZero.csv",header=TRUE,row.names=1)
test<-read.csv(simMatrix.csv,header=TRUE,row.names=1)
#test<-read.csv("fake.csv",header=TRUE,row.names=1)
testMat<-as.matrix(test)



#CHOOSE THE OUTPUT FILE FORMAT
png(paste(outDir,"/conservation_map.png",sep=""), width=4000, height=4000, antialias="default")  #if you have too many features you might have withe stripes in the figure. to avoid this either choose the pdf as output,
#bmp(paste(outDir,"/conservation_map.bmp",sep=""), width=850, height=850, antialias="default") 
#pdf("heatMap2.pdf")                       #either play with the graphical parameters (width height res....usually enlarging the figure should fix the problem..or just use the default)



#CHOOSE THE COLOR (custom, default or RColorBrew), THE SCALE, THE COLUMN ORDER
library(gplots)                    #allows to use heatmap.2 function   
#library("ALL")
#load list of colors. Use it if you do not want to write your own array of colors. This library contines list of colors with all the nuances. That is are scaled. The names are:
#Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd 
library(RColorBrewer) 
#mypalette<-brewer.pal(9,"YlGnBu") #color brewer example: this means that you create an array called mypalette having 9 colors taken from the YlGnBu color scale
#colTest=c("Slate Blue","pink")    #example of custom color array
#data("ALL")                       #load some data...examples...is useless


#HEATMAP DOES NOT WORK IF THERE IS A SINGLE GENOME (A COLUMN). THEREFORE THE FUNCTION IMAGE NEEDS TO BE USED INSTEAD
if (ncol(testMat) == 1) { 
	sortedMatrix=testMat[order(testMat[,1],decreasing=F),]                # sort the matrix in a descending order
	sortedMatrixMat<-as.matrix(sortedMatrix)
	image(t(sortedMatrixMat) , col=colorpanel(40, "red", "black", "green"),xlab=colnames(test),axes=F)            # you need to invert x and y axis, that is swap the row with the columns with the function t()
	axis(2, at=seq(0,1, length.out=4) , labels=rownames(sortedMatrixMat),par(cex=0.5), las=1)            #"at" is telling you where to put the axis. 
                                                                                                             #"seq" indicates that you wanna have intervals. 0,1 means from the bottom left to the top left. 
													     #"length.out=4" means that there will be 4 of them
	quit()
}


#"ward", "single", "complete", "average", "mcquitty", "median" or "centroid"

#OUTPUT THE CLUSTER PRODUCED BY THE HIRARCHICAL CLUSTERING
###WARNING### It is very usefult putting the heatmap values to y or x so that you can see the objects that are inside and check the values

dd.col  <- as.dendrogram(hclust(dist(t(testMat))))
dd.row  <- as.dendrogram(hclust(dist(testMat),method="single"))
dd.col.rev=rev(dd.col)


myCol <- c("black" , "SkyBlue4" , "SkyBlue3" , "SkyBlue2" , "SkyBlue1" , "SkyBlue" , "gold3" , "gold2" , "gold1" , "gold" , "firebrick4" , "firebrick3" ,"firebrick2" ,"firebrick1" )
myBreaks <- c(-1 , 0 , 0.1 , 0.25 , 0.5 , 0.75 , 1 , 1.25 , 1.5 , 1.75 , 2 , 25 , 50 , 75 , 100)
#see bars!
#myCol <- c("black", "SkyBlue4", "SkyBlue3", "SkyBlue2", "SkyBlue1" , "SkyBlue")
#myBreaks <- c(-1 , 0 , 0.1 , 0.25 , 0.5 , 0.75 , 1 )

#par(mar=c(0,0,0,0))

#I suppress the plotting of row dendrogram (although the dendrogram is used to order the rows) because if there are too many leaves it crashes
#x <- heatmap.2(testMat, Colv=dd.col.rev , Rowv=dd.row , col = myCol , breaks = myBreaks , scale="none",key=TRUE,keysize=0.7,symkey=FALSE,density.info="none",trace="none",cexRow=0.2,cexCol=3)
x <- heatmap.2(testMat, Colv=dd.col.rev , Rowv=dd.row , col = myCol , breaks = myBreaks , scale="none",key=FALSE, keysize=0.5 , density.info="none",trace="none", cexRow=0.2,cexCol=6 , margins = c(70, 70))







dev.off()









#revC = identical(Colv, "Rowv"),
#revC=TRUE
#indColumns <- rowInd(testMat, na.rm=TRUE)
#print(indColumns)



#integ=c(2,1)
#heatmap(testMat,scale="row",col=colTest ,Colv=integ)
#heatmap(testMat,scale="row",col=colTest)

#write.table(hc2Newick(x$rowDendrogram), file="/users/cn/gbussotti/Desktop/pig2/results/ibs_new.txt",row.names=FALSE,col.names=FALSE)
