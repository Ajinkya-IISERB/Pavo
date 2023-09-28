#### Region 1

colvec <- c("tomato4","darkolivegreen","turquoise3","peru","orange2","darkseagreen1","slateblue4","cornflowerblue","burlywood4","brown4",
"blueviolet","bisque4","aquamarine1","aquamarine3","aquamarine4","azure3","azure4","brown3","burlywood3","cadetblue2","chartreuse1","chocolate2","coral","coral2","coral3","coral4","cornsilk4","cyan3","cyan4","darkblue","darkgoldenrod","darkgoldenrod2","darkgoldenrod4","darkmagenta","yellow4","tan","tan1","tan2","tan3","tomato",   "turquoise4", "violetred", "violetred4", "wheat3","plum1","tomato3","hotpink", "wheat4","slategrey", "slateblue","slategray4","slateblue2","slategray","yellowgreen","chartreuse4","cadetblue4","violet","navy","khaki","gray50","slategray1","turquoise","olivedrab1","slateblue3","violetred2","slategray2","violetred3","slategray3")
callcol <- c("green","skyblue", "red", "yellow", "magenta")

call_pc <- read.table("../PC_PMN1/PC_PMN1_callable_status.bed")
call_pm1 <- read.table("../PM_PMN1/PM_PMN1_callable_status.bed")
het_pm <- read.table("../PM_PMN1/PM_PMN1_het_5k.bed")
het_pc <- read.table("../PC_PMN1/PC_PMN1_het_5k.bed")

scaf <- "chr29"
pdf(paste(scaf,"_int_dist.pdf",sep=""),width=12,height=9)
par(mar=c(4,6,2.5,2.5))
read.table(file="PM_PMN1.un.dc.txt",header=FALSE)->M
M[M$V1==scaf,]-> Ms
#subset(M, V1='chr29' & V2>35375600 & V3<35917800)->Ms
plot(c(30000000:31000000),rep(5,(31000000-30000000)+1), type="n",ylim=c(0,20),main=scaf,xlab="",ylab="",axes=FALSE,cex.lab=3,col.lab="red")
axis(1, xlim=c(0,max(Ms$V3)+1),col="black",las=1,lwd=4,font.axis=2)
mtext(side=1, line=2, "Genomic Co-ordinates", col="black", font=2,cex=1.25)
call_pc[call_pc$V1==scaf,] -> Cs
rect(Cs$V2,rep(0,length(Cs$V2)),Cs$V3,rep(2,length(Cs$V2)),col=callcol[as.factor(Cs$V4)])
call_pm1[call_pm1$V1==scaf,] -> Cs
rect(Cs$V2,rep(2,length(Cs$V2)),Cs$V3,rep(4,length(Cs$V2)),col=callcol[as.factor(Cs$V4)])
rect(Ms$V2,rep(4,length(Ms$V2)),Ms$V3,rep(6,length(Ms$V2)),col=colvec[Ms$V4])
read.table(file="PC_PMN1.un.dc.txt",header=FALSE)->M
M[M$V1==scaf,]-> Ms
rect(Ms$V2,rep(6.1,length(Ms$V2)),Ms$V3,rep(8,length(Ms$V2)),col=colvec[Ms$V4])

het_pm[het_pm$V1==scaf,] -> Hs
Hs$V7 <- log(Hs$V4+1)
Hs$V8 <- Hs$V7+8.1-min(Hs$V7)
lines(Hs$V3,Hs$V8,type="l",col="blue")
het_pc[het_pc$V1==scaf,] -> Hs
Hs$V7 <- log(Hs$V4+1)
Hs$V8 <- Hs$V7+10.1-min(Hs$V7)
lines(Hs$V3,Hs$V8,type="l",col="blue")

gene <- read.table("../psuedo/tmp.bed")
rect(gene$V2,rep(12,length(gene$V2)),gene$V3,rep(14,length(gene$V2)),col="blue")
text(30792204,14.5,labels="ADAMTSL1")
rect(30792204,rep(12.8,length(gene$V2)),30968165,rep(13,length(gene$V2)),col="blue")
text(30446370,14.5,labels="SH3GL2")
rect(30446370,rep(12.8,length(gene$V2)),30468523,rep(13,length(gene$V2)),col="blue")
text(30179252,14.5,labels="CNTLN")
rect(30179252,rep(12.8,length(gene$V2)),30347714,rep(13,length(gene$V2)),col="blue")


legend("top",legend=seq(0,83,1),ncol = 12,fill=colvec,cex=0.7,bty="n",text.font=2)
text(max(Ms$V3)/2,21,"Atomic Intervals",cex=1,font=2)
legend("top",inset=0.2,ncol = 3,legend=c(unique(Cs$V4)), fill=callcol,cex=0.7,bty="n",text.font=2)
#legend("topleft",inset=.05,fill="blue", cex=1.5, bty="n", legend=expression(bold(paste("Heterozygosity ", theta[w]))))
y.label.position<-c(1,3,5,7,9,11,13)
y.label <- c("PC call", "PM Call", "PM AI", "PC AI", "PM het", "PC het","Genes")
axis(2, at=y.label.position, label=y.label,las=2, lwd=0,font=2,cex.axis=1.35)
box(lwd=2)
dev.off()
}

####################################################################################################
#### Region 2

colvec <- c("tomato4","darkolivegreen","turquoise3","peru","orange2","darkseagreen1","slateblue4","cornflowerblue","burlywood4","brown4",
"blueviolet","bisque4","aquamarine1","aquamarine3","aquamarine4","azure3","azure4","brown3","burlywood3","cadetblue2","chartreuse1","chocolate2","coral","coral2","coral3","coral4","cornsilk4","cyan3","cyan4","darkblue","darkgoldenrod","darkgoldenrod2","darkgoldenrod4","darkmagenta","yellow4","tan","tan1","tan2","tan3","tomato",   "turquoise4", "violetred", "violetred4", "wheat3","plum1","tomato3","hotpink", "wheat4","slategrey", "slateblue","slategray4","slateblue2","slategray","yellowgreen","chartreuse4","cadetblue4","violet","navy","khaki","gray50","slategray1","turquoise","olivedrab1","slateblue3","violetred2","slategray2","violetred3","slategray3")
callcol <- c("green","skyblue", "red", "yellow", "magenta")

call_pc <- read.table("../PC_PMN1/PC_PMN1_callable_status.bed")
call_pm1 <- read.table("../PM_PMN1/PM_PMN1_callable_status.bed")
het_pm <- read.table("../PM_PMN1/PM_PMN1_het_5k.bed")
het_pc <- read.table("../PC_PMN1/PC_PMN1_het_5k.bed")

scaf <- "chr26"
pdf(paste(scaf,"_int_dist.pdf",sep=""),width=12,height=9)
par(mar=c(4,6,2.5,2.5))
read.table(file="PM_PMN1.un.dc.txt",header=FALSE)->M
M[M$V1==scaf,]-> Ms
#subset(M, V1='chr29' & V2>35375600 & V3<35917800)->Ms
plot(c(57891100:58041100),rep(5,(58041100-57891100)+1), type="n",ylim=c(0,20),xlab="",ylab="",axes=FALSE,cex.lab=3,col.lab="red")
axis(1, xlim=c(0,max(Ms$V3)+1),col="black",las=1,lwd=4,font.axis=2)
mtext(side=1, line=2, "Genomic Co-ordinates", col="black", font=2,cex=1.25)

call_pc[call_pc$V1==scaf,] -> Cs
rect(Cs$V2,rep(0,length(Cs$V2)),Cs$V3,rep(2,length(Cs$V2)),col=callcol[as.factor(Cs$V4)])
call_pm1[call_pm1$V1==scaf,] -> Cs
rect(Cs$V2,rep(2,length(Cs$V2)),Cs$V3,rep(4,length(Cs$V2)),col=callcol[as.factor(Cs$V4)])
rect(Ms$V2,rep(4,length(Ms$V2)),Ms$V3,rep(6,length(Ms$V2)),col=colvec[Ms$V4])
read.table(file="PC_PMN1.un.dc.txt",header=FALSE)->M
M[M$V1==scaf,]-> Ms
rect(Ms$V2,rep(6.1,length(Ms$V2)),Ms$V3,rep(8,length(Ms$V2)),col=colvec[Ms$V4])
#abline(v=57958756,lty=2,lwd=3)
#abline(v=57988172,lty=2,lwd=3)

het_pm[het_pm$V1==scaf,] -> Hs
Hs$V7 <- log(Hs$V4+1)
Hs$V8 <- Hs$V7+7.1-min(Hs$V7)
lines(Hs$V3,Hs$V8,type="l",col="green")
het_pc[het_pc$V1==scaf,] -> Hs
Hs$V7 <- log(Hs$V4+1)
Hs$V8 <- Hs$V7+8.1-min(Hs$V7)
lines(Hs$V3,Hs$V8,type="l",col="blue")

gene <- read.table("../psuedo/reg2_exons.bed")
rect(gene$V2,rep(12,length(gene$V2)),gene$V3,rep(14,length(gene$V2)),col="blue")
rect(min(gene$V2),rep(12.8,length(gene$V2)),max(gene$V3),rep(13,length(gene$V2)),col="blue")
text(min(gene$V2),14.5,labels="CALD1")

legend("top",legend=seq(0,83,1),ncol = 12,fill=colvec,cex=0.7,bty="n",text.font=2)
text(57965000,17,"Atomic Intervals",cex=1,font=2)
legend("top",inset=0.2,ncol = 3,legend=c(unique(Cs$V4)), fill=callcol,cex=0.7,bty="n",text.font=2)
#legend("topleft",inset=.05,fill="blue", cex=1.5, bty="n", legend=expression(bold(paste("Heterozygosity ", theta[w]))))
y.label.position<-c(1,3,5,7,9,11,13)
y.label <- c("PC call", "PM Call", "PM AI", "PC AI", "PM het", "PC het","Genes")
axis(2, at=y.label.position, label=y.label,las=2, lwd=0,font=2,cex.axis=1.35)
box(lwd=2)
dev.off()
}
