options(expressions=500000) #this option is needed if you have many features in your matrix and you have to tell R that it is not an infinite loop, but just a bit longer
###example commandLine: R --slave --args simMatrix.csv experiments/exp_1/results/ < qualityCheck/utility/heatmap.R
Args <- commandArgs();
simMatrix.csv=Args[4]
outDir=Args[5]


#LIBRARY NEEDED IF THE HEATMAP IS HUGE (i.e. 17000 transcripts)
# if the package is not installed use the command  install.packages("fastcluster") to install it
library('fastcluster')
library(gplots)                    #allows to use heatmap.2 function   
library(RColorBrewer) 



#CHOOSE THE INPUT FILE
#test<-read.csv("heatMap_noAllZero.csv",header=TRUE,row.names=1)
test<-read.csv(simMatrix.csv,header=TRUE,row.names=1)
#test<-read.csv("fake.csv",header=TRUE,row.names=1)

#banded
for (y in 1:length(colnames(test))) {for (x in 1:length(test[,y])) { if (test[x,y] > 10){test[x,y] = 10}   }  }
#for (y in 1:length(colnames(test))) {for (x in 1:length(test[,y])) { if (test[x,y] == 0){}else{test[x,y] = log(test[x,y])}   }  }

#remove rows all 0
toRemove=c()
for (x in 1:length(rownames(test))) { counter=0 ;  for (y in 1:length(colnames(test))) { if (test[x,y] == 0){counter=counter+1 }}; if (counter==6){toRemove=c(toRemove,x)}  }
test=test[-toRemove,]


#to matrix
testMat<-as.matrix(test)




#CHOOSE THE OUTPUT FILE FORMAT
png(paste(outDir,"/conservation_map.png",sep=""), width=1000, height=1000, antialias="default")  #if you have too many features you might have withe stripes in the figure. to avoid this either choose the pdf as output,


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

#OUTPUT THE CLUSTER PRODUCED BY THE HIRARCHICAL CLUSTERING
dd.col  <- as.dendrogram(hclust(dist(t(testMat))))
dd.row  <- as.dendrogram(hclust(dist(testMat) , method="complete"))
dd.col.rev=rev(dd.col)




#I suppress the plotting of row dendrogram (although the dendrogram is used to order the rows) because if there are too many leaves it crashes
#x <- heatmap.2(testMat, Colv=dd.col.rev , Rowv=dd.row ,col=colorpanel(40, "blue", "yellow", "green"),scale="none",key=TRUE,keysize=0.7,symkey=FALSE,density.info="none",trace="none",cexRow=0.5,cexCol=1.1 , dendrogram="column")


myCol <- c("blue", "green", "yellow", "orange", "red")
myBreaks <- c(0 , 0.5 , 1 , 5 , 50 , 10000000)
#x <- heatmap.2(testMat, Colv=dd.col.rev , Rowv=dd.row , col = myCol , breaks = myBreaks , trace="none",cexRow=0.5,cexCol=1.1 , dendrogram="column")
#legend("left", fill = myCol, legend = c("1", "2 to 5", "6 to 50", ">50"))

x <- heatmap.2(testMat, Colv=dd.col.rev , Rowv=dd.row  , trace="none",cexRow=0.5,cexCol=1.1 , dendrogram="column")
dev.off()












#OTHER CMDS
#library("ALL")
#load list of colors. Use it if you do not want to write your own array of colors. This library contines list of colors with all the nuances. That is are scaled. The names are:
#Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd 
#mypalette<-brewer.pal(9,"YlGnBu") #color brewer example: this means that you create an array called mypalette having 9 colors taken from the YlGnBu color scale
#colTest=c("Slate Blue","pink")    #example of custom color array
#data("ALL")                       #load some data...examples...is useless



###WARNING### It is very usefult putting the heatmap values to y or x so that you can see the objects that are inside and check the values
#y <- heatmap.2(testMat, scale="none")
#dend.col <- y$colDendrogram         #dend.col is an object of the class dendrogram
#rev.dend.col=rev(dend.col)          #use rev() to reverse a dendrogram. Then I recompute another heatmap imposing for the column the reversed dendrogram (option Colv)

#revC = identical(Colv, "Rowv"),
#revC=TRUE
#indColumns <- rowInd(testMat, na.rm=TRUE)
#print(indColumns)

#pdf("heatMap2.pdf")   #either play with the graphical parameters (width height res....usually enlarging the figure should fix the problem..or just use the default)


#integ=c(2,1)
#heatmap(testMat,scale="row",col=colTest ,Colv=integ)
#heatmap(testMat,scale="row",col=colTest)

#write.table(hc2Newick(x$rowDendrogram), file="/users/cn/gbussotti/Desktop/pig2/results/ibs_new.txt",row.names=FALSE,col.names=FALSE)


#hc <- hclust(dist(testMat), method="complete")   #you can make hirarchical clusteing of a matrix just doing this. this is what the heatmap function does to compute the dendrogram
#library(ctc)                                      #library useful for hanling clustering

#hclust2treeview(testMat,file="cluster.cdt")    #treeview readable output
#newick_out<-hc2Newick(hc, flat=TRUE)            #newick format output
#print(newick_out)


#PRINT SOME USEFUL NUMBERS  (only if you use "scale", as done by default)
#mean_rows <- rowMeans(testMat, na.rm=TRUE)
#print(mean_rows)
#mean_col <- colMeans(testMat, na.rm=TRUE)
#print(mean_col)
