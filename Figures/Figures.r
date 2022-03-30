pc <- read.table("PC_PC1/PC_PC1.un.0.txt")
ps <- read.table("psuedo/PM_PMN1.PC_PMN1.merged.0.txt")
pm <- read.table("PM_PMN1/PM_PMN1.un.0.txt")
pc_all = list.files("PC_PC1",pattern="PC_PC1.combined.*.txt", full.names = TRUE)
ps_all = list.files("psuedo",pattern="PM_PMN1.PC_PMN1.merged.combined.*.txt", full.names = TRUE)
pm_all = list.files("PM_PMN1",pattern="PM_PMN1.combined.*.txt", full.names = TRUE)
d <- read.table("atm_sea_data.txt",header=T)

jpeg(file="Pavo_psmc_with_temp_bs.jpeg",width = 14000, height = 5000,res=300)
par(mar=c(20.1,23.1,8.1,22.1),mgp = c(18, 10, 0))
xend <- 6.6
plot(log10(ps$V1),ps$V2,type="s",col="red",pch=16,xaxt="n",yaxt="n",xlab="Time (MYA)",ylab=substitute(paste("Effective population size ",italic('N'),"e (x 10"^"4",")",sep="")),cex.axis=5,lwd=5,cex.lab=5,bty="n",xlim=c(4,xend),ylim=c(0,20)
,font=2,font.lab=2,font.axis=5)

axis(2,at=c(0,5,10,15,20),labels=c(0,5,10,15,20),cex.axis=5,lwd=10,font=2)
axis(1,at=log10(c(1,10000,20000,30000,40000,50000,100000,200000,400000,600000,800000,1000000,
1500000,2000000,3000000,4000000,5000000)),labels=c(0,0.01,0.02,0.03,0.04,0.05,0.1,0.2,
0.4,0.6,0.8,1,1.5,2,3,4,5),cex.axis=5,lwd=10,font=2)

x<-10000
Csze<-5
den=3
rect(log10(127000),170000/x,log10(1), 20, col=scales::alpha(rgb(col2rgb("cyan")[1,],col2rgb("cyan")[2,],col2rgb("cyan")[3,],max = 255), 0.25),border = "transparent")
text(log10(3000),180000/x, "Holocene",cex=Csze,font=2)
rect(log10(110000),170000/x,log10(11700), 20, col=scales::alpha(rgb(col2rgb("black")[1,],col2rgb("black")[2,],col2rgb("black")[3,],max = 255), 0.25),border = "transparent")
text(log10(35599),180000/x, "Last glacial period",cex=Csze,font=2)

rect(log10(2580000),170000/x,log10(128000),20,  col=scales::alpha(rgb(col2rgb("yellow")[1,],col2rgb("yellow")[2,],col2rgb("yellow")[3,],max = 255), 0.25),border = "transparent")
text(log10(700000),180000/x, "Pleistocene",cex=Csze,font=2)

rect(log10(5330000),170000/x,log10(2580000),20,  col=scales::alpha(rgb(col2rgb("green")[1,],col2rgb("green")[2,],col2rgb("green")[3,],max = 255), 0.25),border = "transparent")
text(log10(3500000),180000/x, "Pleiocene",cex=Csze,font=2)

for (i in pc_all){
read.table(file=i,header=FALSE)->a
lines(log10(a$V1),a$V2, type="s",col=scales::alpha(rgb(col2rgb("blue")[1,],col2rgb("blue")[2,],col2rgb("blue")[3,],max = 255), 0.25))
}
for (i in pm_all){
read.table(file=i,header=FALSE)->a
lines(log10(a$V1),a$V2, type="s",col=scales::alpha(rgb(col2rgb("darkgreen")[1,],col2rgb("darkgreen")[2,],col2rgb("darkgreen")[3,],max = 255), 0.25))
}
for (i in ps_all){
read.table(file=i,header=FALSE)->a
lines(log10(a$V1),a$V2, type="s",col=scales::alpha(rgb(col2rgb("red")[1,],col2rgb("red")[2,],col2rgb("red")[3,],max = 255), 0.25))
}

points(log10(pc$V1),pc$V2,type="s",col="blue",lwd=15)
points(log10(pm$V1),pm$V2,type="s",col="darkgreen",lwd=15)
points(log10(ps$V1),ps$V2,type="s",col="red",lwd=15)

par(new = T,mar=c(25.1,25.1,12.1,27.1),mgp = c(40, 12, 3))
myPlot <- function(x,y,zcol) {
    plot(x, y, yaxt = "n",xaxt="n",type="l",col=zcol,xlim=c(4,xend),xlab="",ylab="",bty="n",lwd=5)
    axis(4,cex.axis=5,lwd=10,font=2)
    ylab=expression(Temperature~(degree*C))
    mtext(ylab, side=4, line = 24,cex=5,font=2)
}

myPlot(log10(d$Time*1000),d$Tsurf,adjustcolor("brown", alpha.f = 0.5))

dev.off()

#######################################################################################

jpeg("Figure_3c.jpeg",res=300,height=2000,width=3000)
par(mar = c(2, 3, 1, 1))

a <- read.table("split-time-est.txt",header=T,sep="\t")

plot(a$pos,a$Split.time.estimate,xaxt="n",xlab="",ylab="",type="n",axes=F,cex.main=1.75,xlim=c(0,0.5))

rect(-1,2.58,12,0, col=scales::alpha(rgb(col2rgb("yellow")[1,],col2rgb("yellow")[2,],col2rgb("yellow")[3,],max = 255), 0.25),border="transparent")
text(0.25,1.7,"Pleistocene",cex=1.5,font=4,col="red")

rect(-1,5.33,12,2.58, col=scales::alpha(rgb(col2rgb("green")[1,],col2rgb("green")[2,],col2rgb("green")[3,],max = 255), 0.25),border="transparent")
text(0.25,3.7,"Pliocene",cex=1.5,font=4,col="red")

rect(-1,7.25,12,5.33, col=scales::alpha(rgb(col2rgb("orange")[1,],col2rgb("orange")[2,],col2rgb("orange")[3,],max = 255), 0.25),border="transparent")
text(0.25,5.7,"Miocene",cex=1.5,font=4,col="red")

axis(2, col="black",las=1,font=2,font.axis=4,cex.axis=1.25)
mtext(side=2, line=1.6, "Time (MYA)", col="black", font=2,cex=1.5)
mtext(side=1, line=0.7, "Split time estimates",col="black", font=2,cex=1.5)

points(a$pos,a$Split.time.estimate,pch=19)
text(0.2,4.2,"CYTB-Kimball et al. (1997)",cex=1.5,font=2)
text(0.2,6.2,"D-loop-Kimball et al. (1997)",cex=1.5,font=2)
text(0.12,1.5,"NM-Jetz et al. (2012)",cex=1.5,font=2)
text(0.28,2.27,"NM-Stein et al. (2015)",cex=1.5,font=2)
text(0.3,2.96,"CYTB-Ouyang et al. (2009)",cex=1.5,font=2)
text(0.11,3.3,"UCE/WGS-Kimball et al. (2019)",cex=1.5,font=2)
text(0.28,2.58,"NM-Wang et al. (2016)",cex=1.5,font=2)
text(0.32,3.26,"UCE-Chen et al. (2021)",cex=1.5,font=2)
text(0.23,2.02,"UCE - This study",cex=1.5,font=2)
rect(0.4,1.1,0.44,3, col=scales::alpha(rgb(col2rgb("red")[1,],col2rgb("red")[2,],col2rgb("orange")[3,],max = 255), 0.5),border="transparent")
text(0.4,1,"PSMC psuedo-diploid - this study",cex=1.4,font=2)
rect(0.4,2.99,0.44,3,col="black")
rect(0.4,1.79,0.44,1.8,col="black")
text(0.48,1.2,"Final split",cex=1.4,font=2)
text(0.48,1.4,"(~ 1.1 MYA)",cex=1.4,font=2)
text(0.48,1.9,"Gene flow",cex=1.4,font=2)
text(0.48,2.1,"(~ 1.8 MYA)",cex=1.4,font=2)
text(0.48,2.9,"Initial split",cex=1.4,font=2)
text(0.48,3.1,"(~ 3 MYA)",cex=1.4,font=2)
box()

dev.off()

library(ggplot2)
library(ggrepel)
a <- read.table("all_pathways_fst0.9.bed")
a[which(a$V9>0.01),]->b
pdf("pathways_fst_boxplots.pdf",width=13,height=10)
ggplot(a,aes(x=V12,y=V9))+ geom_boxplot(outlier.shape=NA,col="red")+geom_jitter(data=b,aes(x=V12,y=V9),position=position_jitter(width=.1, height=0),col="blue")+ggrepel::geom_text_repel(data = b, aes(label = V11))+labs(x="Signalling Pathways", y="Fixed site density")+coord_cartesian(ylim = c(0,0.030))
dev.off()

pdf("pathways_fst_boxplots.pdf",width=13,height=10)
ggplot(a,aes(x=V12,y=V9))+ geom_boxplot(outlier.shape=NA,col="red")+geom_jitter(data=b,aes(x=V12,y=V9),position=position_jitter(width=.1, height=0),col="blue")+ggrepel::geom_text_repel(data = b, aes(label = V11))+labs(x="Signalling Pathways", y="Fixed site density")+coord_cartesian(ylim = c(0,0.030))+
scale_x_discrete(labels=c("ACTIN CYTO.","ECM INT.","FOCAL ADHESION","HEDGEHOG","MELANOGENESIS","NOTCH","TGF−BETA","WNT"))

jpeg("pathways_fst_boxplots.jpeg",width=13,height=10,units="in",res=1000)
ggplot(a,aes(x=V12,y=V9))+ geom_boxplot(outlier.shape=NA,col="red")+geom_jitter(data=b,aes(x=V12,y=V9),position=position_jitter(width=.1, height=0),col="blue")+ggrepel::geom_text_repel(data = b, aes(label = V11),size = 5)+labs(x="Signalling Pathways", y="Fixed site density")+coord_cartesian(ylim = c(0,0.030))+
scale_x_discrete(labels=c("ACTIN CYTO.","ECM INT.","FOCAL ADHESION","HEDGEHOG","MELANOGENESIS","NOTCH","TGF-BETA","WNT"))+theme(axis.text=element_text(size=14,color="black"),axis.title=element_text(size=16,face="bold"))
dev.off()


