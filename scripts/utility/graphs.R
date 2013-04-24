#R histogram plot
png("mouseMapping.png", width=1000, height=900 , res=120)
mouse=read.table("mouse_mapping_count",header=FALSE)
par(ps=10)
mouse_hist=hist(mouse$V2,  plot=F , breaks=seq(-1,max(mouse$V2),by=1))
clr <- rep("orange", length(mouse_hist$counts));
clr[1] <- "black";
plot(mouse_hist,  col=clr, xlab="bin",ylab="count", main="",xlim=c(-1,20),axes=F, label=T)
axis(1,at=mouse_hist$mids,labels=0:max(mouse$V2))
axis(2)
dev.off()

png("humanMapping.png", width=1000, height=900 , res=120)
human=read.table("human_mapping_count",header=FALSE)
par(ps=10)

human_hist=hist(human$V2, plot = F , breaks=seq(-1,max(human$V2),by=1))
clr <- rep("red", length(human_hist$counts));
clr[1] <- "black";
plot(human_hist,  col=clr, xlab="number of homologs",ylab="count", main="",xlim=c(-1,20),axes=F, label=T)
axis(1,at=human_hist$mids,labels=0:max(human$V2))
axis(2)
dev.off()

print("mouse summary")
summary(mouse$V2)
print("human summary")
summary(human$V2)






#repeats
png("repeatCoverage.png", width=1000, height=900 , res=120)

m=read.table("mouse_repeatCoverage")
mdat<-m$V2
density_mdat<-density(mdat)
plot(density_mdat,xlim=c(0,100),xaxs="i",col="orange",xlab="fraction covered by repeats",main="")

h=read.table("human_repeatCoverage")
hdat<-h$V2
density_hdat<-density(hdat)
lines(density_hdat,xlim=c(0,100),xaxs="i",col="red")
legend("topright" , c("human","mouse") , fill=c("red","orange"))

dev.off()





#exon Conservation
png("exonConservation.png", width=1500, height=900 , res=120)
exonConservation=read.table("exonNumberConservation")
par(ps=10)
min_br=min(exonConservation$V5) -1 
exConsHist=hist(exonConservation$V5, plot = F , breaks=seq(min_br,max(exonConservation$V5),by=1))
clr <- rep("gray", length(exConsHist$counts));
clr[16] <- "green";
plot(exConsHist,  col=clr,xlab="bin",ylab="count",axes=F, main="" , label=T )
axis(1,at=exConsHist$mids, labels=seq(-15 ,13))
axis(2)
dev.off()










