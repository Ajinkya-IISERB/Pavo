#Make single base pair window over genome
bedtools makewindows -g ../Pmut.genome -w 1 -s 1 > Pmut.1bp.bed
#Make 1Kb windows over genome
bedtools makewindows -g ../Pmut.genome -w 1000 -s 1000 > Pmut.1Kbp.bed

#get genomic regions for each individual PSMC file.
#P.cristatus
decode2bed.pl -s100 PC_PMN1.un.psmc > PC_PMN1.un.dc.txt
#P. muticus individual 1
decode2bed.pl -s100 PM_PMN1.un.psmc > PM_PMN1.un.dc.txt
#P. muticus individual 2
decode2bed.pl -s100 PM_PMN2.un.psmc > PM_PMN2.un.dc.txt
#P. muticus individual 3
decode2bed.pl -s100 PM_PMN3.un.psmc > PM_PMN3.un.dc.txt
#P. muticus individual 4
decode2bed.pl -s100 PM_PMN4.un.psmc > PM_PMN4.un.dc.txt

#P. cristatus and P. muticus
#get intersections between all the atomic intervals between pai of individuals.
bedtools intersect -a PM_PMN1.un.dc.txt -b PC_PMN1.un.dc.txt -wao | awk '$6!="." {print $0}' | awk '{print $6,$7,$8,$4,$9}' | sed 's/ /\t/g' > PM_PC_DC.txt
#Map onto 1 bp and 1KB windows to get DAI scores.
bedtools intersect -a Pmut.1bp.bed -b PM_PC_DC.txt -wao | awk '{print $1,$2,$3,$7,$8,$8-$7}' | sed 's/ /\t/g' | bedtools map -a ../Pmut.1Kbp.bed -b stdin -c 6 -o sum > PM_PC_DAI.bed
#P. muticus ind1 and P. muticus ind2
#get intersections between all the atomic intervals between pai of individuals.
bedtools intersect -a PM_PMN1.un.dc.txt -b PC_PMN2.un.dc.txt -wao | awk '$6!="." {print $0}' | awk '{print $6,$7,$8,$4,$9}' | sed 's/ /\t/g' > PM1_PM2_DC.txt
#Map onto 1 bp and 1KB windows to get DAI scores.
bedtools intersect -a Pmut.1bp.bed -b PM1_PM2_DC.txt -wao | awk '{print $1,$2,$3,$7,$8,$8-$7}' | sed 's/ /\t/g' | bedtools map -a ../Pmut.1Kbp.bed -b stdin -c 6 -o sum > PM1_PM2_DAI.bed

for i in PM_PMN1 PM_PMN2 PM_PMN3 PM_PMN4 PM_PMN5 PM_PMN6
do
for y in PM_PMN1 PM_PMN2 PM_PMN3 PM_PMN4 PM_PMN5 PM_PMN6
do
bedtools intersect -a "$i".un.dc.txt -b "$y".un.dc.txt -wao | awk '$6!="." {print $0}' | awk '{print $6,$7,$8,$4,$9}' | sed 's/ /\t/g' > "$i"_"$y"_DC.txt
done
done

for i in PM_PMN1_PM_PMN3_DC.txt PM_PMN1_PM_PMN4_DC.txt PM_PMN1_PM_PMN5_DC.txt PM_PMN1_PM_PMN6_DC.txt
do
y=`echo $i | sed 's/_DC.txt//g'`
bedtools intersect -a Pmut.1bp.bed -b $i -wao | awk '{print $1,$2,$3,$7,$8,$8-$7}' | sed 's/ /\t/g' | bedtools map -a ../Pmut.1Kbp.bed -b stdin -c 6 -o sum > "$y"_DAI.bed
done

############################################
#In R
##############################
#DAI plots
read.table(file="PM_PC_DAI.bed",header=FALSE,fill=TRUE)->M

pdf("PM_PC_DAI_all_chr.pdf")
par(mfrow=c(2,2))
for(i in unique(M$V1))
{
plot(M$V2[M$V1==i],M$V4[M$V1==i],pch=20,xlab="position",main=i,ylab="DAI")
}
dev.off()
