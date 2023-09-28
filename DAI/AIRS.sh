#get genomic regions for each individual PSMC file.
#P.cristatus
mod_dec.pl PC_PMN1.un.psmc | awk '{ print >$4".pc"}'
#P. muticus individual 1
mod_dec.pl PM_PMN1.un.psmc | awk '{ print >$4".pm1"}'
#P. muticus individual 2
mod_dec.pl PM_PMN2.un.psmc | awk '{ print >$4".pm2"}'
#P. muticus individual 3
mod_dec.pl PM_PMN3.un.psmc | awk '{ print >$4".pm3"}'
#P. muticus individual 4
mod_dec.pl PM_PMN4.un.psmc | awk '{ print >$4".pm4"}'

#Get intersections between P. muticus individual-1 and P. cristatus.
for j in *.pm1
do 
a=`echo "$j" | sed 's/.pm1//g'` 
for i in *.pc
do
y=`echo $i | sed 's/.pc//g'`
z=`bedtools intersect -a "$j" -b "$i" | awk '{s+=$3-$2}END{print s}'`
echo $y $z | sed 's/ /\t/g' >> "$j"_pc_int.txt
done
done
#Get lengths of the intersections.
for i in *.pm1_pc_int.txt
do
y=`echo $i | sed 's/.pm1_pc_int.txt//g'`
cat $i | awk -v a="$y" '$2!="" {print a"\t"$0}' >> pm1_pc_ai_distr.txt
done

########################################
#Get intersections between P. muticus individual-1 and individual-2
for j in *.pm1
do 
a=`echo "$j" | sed 's/.pm1//g'` 
for i in *.pm2
do
y=`echo $i | sed 's/.pm2//g'`
z=`bedtools intersect -a "$j" -b "$i" | awk '{s+=$3-$2}END{print s}'`
echo $y $z | sed 's/ /\t/g' >> "$j"_pm2_int.txt
done
done
#Get lengths of the intersections.
for i in *.pm1_pm2_int.txt
do
y=`echo $i | sed 's/.pm1_pm2_int.txt//g'`
cat $i | awk -v a="$y" '$2!="" {print a"\t"$0}' >> pm1_pm2_ai_distr.txt
done

####################################
#Get intersections between P. muticus individual-1 and individual-3
for j in *.pm1
do 
a=`echo "$j" | sed 's/.pm1//g'` 
for i in *.pm3
do
y=`echo $i | sed 's/.pm3//g'`
z=`bedtools intersect -a "$j" -b "$i" | awk '{s+=$3-$2}END{print s}'`
echo $y $z | sed 's/ /\t/g' >> "$j"_pm3_int.txt
done
done
#Get lengths of the intersections.
for i in *.pm1_pm3_int.txt
do
y=`echo $i | sed 's/.pm1_pm3_int.txt//g'`
cat $i | awk -v a="$y" '$2!="" {print a"\t"$0}' >> pm1_pm3_ai_distr.txt
done

####################################
#Get intersections between P. muticus individual-1 and individual-4
for j in *.pm1
do 
a=`echo "$j" | sed 's/.pm1//g'` 
for i in *.pm4
do
y=`echo $i | sed 's/.pm4//g'`
z=`bedtools intersect -a "$j" -b "$i" | awk '{s+=$3-$2}END{print s}'`
echo $y $z | sed 's/ /\t/g' >> "$j"_pm4_int.txt
done
done
#Get lengths of the intersections.
for i in *.pm1_pm4_int.txt
do
y=`echo $i | sed 's/.pm1_pm4_int.txt//g'`
cat $i | awk -v a="$y" '$2!="" {print a"\t"$0}' >> pm1_pm4_ai_distr.txt
done
########################################################

#Get intersections between P. muticus individual-1 and P. cristatus
for j in *.pm1
do 
a=`echo "$j" | sed 's/.pm1//g'` 
for i in *.pc
do
y=`echo $i | sed 's/.pc//g'`
z=`bedtools intersect -a "$j" -b "$i" | awk '{s+=$3-$2}END{print s}'`
echo $y $z | sed 's/ /\t/g' >> "$j"_pc_int.txt
done
done
#Get lengths of the intersections.
for i in *.pm1_pc_int.txt
do
y=`echo $i | sed 's/.pm1_pc_int.txt//g'`
cat $i | awk -v a="$y" '$2!="" {print a"\t"$0}' >> pm1_pc_ai_distr.txt
done

########################################
#Get intersections between P. muticus individual-2 and P. cristatus
for j in *.pm2
do 
a=`echo "$j" | sed 's/.pm2//g'` 
for i in *.pc
do
y=`echo $i | sed 's/.pc//g'`
z=`bedtools intersect -a "$j" -b "$i" | awk '{s+=$3-$2}END{print s}'`
echo $y $z | sed 's/ /\t/g' >> "$j"_pc_int.txt
done
done
#Get lengths of the intersections.
for i in *.pm2_pc_int.txt
do
y=`echo $i | sed 's/.pm2_pc_int.txt//g'`
cat $i | awk -v a="$y" '$2!="" {print a"\t"$0}' >> pm2_pc_ai_distr.txt
done

####################################
#Get intersections between P. muticus individual-3 and P. cristatus
for j in *.pm3
do 
a=`echo "$j" | sed 's/.pm3//g'` 
for i in *.pc
do
y=`echo $i | sed 's/.pc//g'`
z=`bedtools intersect -a "$j" -b "$i" | awk '{s+=$3-$2}END{print s}'`
echo $y $z | sed 's/ /\t/g' >> "$j"_pc_int.txt
done
done
#Get lengths of the intersections.
for i in *.pm3_pc_int.txt
do
y=`echo $i | sed 's/.pm3_pc_int.txt//g'`
cat $i | awk -v a="$y" '$2!="" {print a"\t"$0}' >> pm3_pc_ai_distr.txt
done

####################################
#Get intersections between P. muticus individual-4 and P. cristatus
for j in *.pm4
do 
a=`echo "$j" | sed 's/.pm4//g'`
for i in *.pc
do
y=`echo $i | sed 's/.pc//g'`
z=`bedtools intersect -a "$j" -b "$i" | awk '{s+=$3-$2}END{print s}'`
echo $y $z | sed 's/ /\t/g' >> "$j"_pc_int.txt
done
done
#Get lengths of the intersections.
for i in *.pm4_pc_int.txt
do
y=`echo $i | sed 's/.pm4_pc_int.txt//g'`
cat $i | awk -v a="$y" '$2!="" {print a"\t"$0}' >> pm4_pc_ai_distr.txt
done

#####################################
#In R
###############################################################
#Figure 1
a <- read.table("pm1_pc_ai_distr.txt")
(a$V3*100)/sum(a$V3)->a$V4
b <- read.table("pm1_pm2_ai_distr.txt")
(b$V3*100)/sum(b$V3)->b$V4

pc <- read.table("PC_PMN1.un.0.txt")
pm1 <- read.table("PM_PMN1.un.0.txt")
diff_ne <- pc$V2-pm1$V2
pm2 <- read.table("PM_PMN2.un.0.txt")
diff_ne_pm <- pm1$V2-pm2$V2

pdf("PC_PM_AI_distr_diff_ne.pdf",width=20,height=10)
layout(matrix(c(1,3,2,4), 2, 2, byrow = TRUE),widths=c(1,1), heights=c(3,1))

par(mar=c(5,5,1,1))
plot(a$V1, a$V2, pch=20, cex=a$V4*5, xlab="Pavo muticus", ylab="Pavo cristatus", font=2, font.lab=4, cex.lab=1.5, cex.main=2, cex.axis=1.5, col=ifelse(a$V4>quantile(a$V4,0.99), "red", ifelse(a$V4>quantile(a$V4,0.95), "blue", ifelse(a$V4>quantile(a$V4,0.9), "green", "black"))))
abline(v=60,lty=2,lwd=2,col="brown")
abline(h=60,lty=2,lwd=2,col="brown")
abline(v=40,lty=2,lwd=2,col="brown")
abline(h=40,lty=2,lwd=2,col="brown")
abline(v=20,lty=2,lwd=2,col="brown")
abline(h=20,lty=2,lwd=2,col="brown")
abline(coef=c(6,1),lty=2,lwd=2,col="red")
abline(coef=c(-6,1),lty=2,lwd=2,col="red")

par(mar=c(5,5,1,1))
plot(diff_ne, type="s", lwd=2, xlab="Atomic Intervals", ylab="Difference in Ne", font=2, font.lab=2, cex.lab=1.4, cex.main=2, cex.axis=1.5, ylim =c(-5,15))
abline(h=0,lty=2,col="red")
abline(h=5,lty=2,col="red")
abline(h=10,lty=2,col="red")

par(mar=c(5,5,1,1))
plot(b$V1, b$V2, pch=20, cex=b$V4*5, xlab="Pavo muticus Indv 1", ylab="Pavo muticus Indv 2", font=2, font.lab=4, cex.lab=1.5, cex.main=2, cex.axis=1.5, col=ifelse(b$V4>quantile(b$V4,0.99), "red", ifelse(b$V4>quantile(b$V4,0.95), "blue", ifelse(b$V4>quantile(b$V4,0.9), "green", "black"))))
abline(v=60,lty=2,lwd=2,col="brown")
abline(h=60,lty=2,lwd=2,col="brown")
abline(v=40,lty=2,lwd=2,col="brown")
abline(h=40,lty=2,lwd=2,col="brown")
abline(v=20,lty=2,lwd=2,col="brown")
abline(h=20,lty=2,lwd=2,col="brown")
abline(coef=c(6,1),lty=2,lwd=2,col="red")
abline(coef=c(-6,1),lty=2,lwd=2,col="red")

par(mar=c(5,5,1,1))
plot(diff_ne_pm, type="s", lwd=2, xlab="Atomic Intervals", ylab="Difference in Ne", font=2, font.lab=2, cex.lab=1.4, cex.main=2, cex.axis=1.5, ylim =c(-5,15))
abline(h=0,lty=2,col="red")
abline(h=5,lty=2,col="red")
abline(h=10,lty=2,col="red")

dev.off()

###############################################################
#Figure S1
a <- read.table("pm1_pc_ai_distr.txt")
(a$V3*100)/sum(a$V3)->a$V4
b <- read.table("pm1_pm2_ai_distr.txt")
(b$V3*100)/sum(b$V3)->b$V4
c <- read.table("pm1_pm3_ai_distr.txt")
(c$V3*100)/sum(c$V3)->c$V4
d <- read.table("pm1_pm4_ai_distr.txt")
(d$V3*100)/sum(d$V3)->d$V4


pdf("4panel_distr_diff_ne.pdf",width=10,height=10)
par(mfrow=c(2,2))

par(mar=c(5,5,1,1))
plot(a$V1, a$V2, pch=20, cex=a$V4*5, xlab="Pavo muticus Indv 1", ylab="Pavo cristatus", font=2, font.lab=4, cex.lab=1.5, cex.main=2, cex.axis=1.5, col=ifelse(a$V4>quantile(a$V4,0.99), "red", ifelse(a$V4>quantile(a$V4,0.95), "blue", ifelse(a$V4>quantile(a$V4,0.9), "green", "black"))))
abline(v=60,lty=2,lwd=2,col="brown")
abline(h=60,lty=2,lwd=2,col="brown")
abline(v=40,lty=2,lwd=2,col="brown")
abline(h=40,lty=2,lwd=2,col="brown")
abline(v=20,lty=2,lwd=2,col="brown")
abline(h=20,lty=2,lwd=2,col="brown")
abline(coef=c(6,1),lty=2,lwd=2,col="red")
abline(coef=c(-6,1),lty=2,lwd=2,col="red")

par(mar=c(5,5,1,1))
plot(b$V1, b$V2, pch=20, cex=b$V4*5, xlab="Pavo muticus Indv 1", ylab="Pavo muticus Indv 2", font=2, font.lab=4, cex.lab=1.5, cex.main=2, cex.axis=1.5, col=ifelse(b$V4>quantile(b$V4,0.99), "red", ifelse(b$V4>quantile(b$V4,0.95), "blue", ifelse(b$V4>quantile(b$V4,0.9), "green", "black"))))
abline(v=60,lty=2,lwd=2,col="brown")
abline(h=60,lty=2,lwd=2,col="brown")
abline(v=40,lty=2,lwd=2,col="brown")
abline(h=40,lty=2,lwd=2,col="brown")
abline(v=20,lty=2,lwd=2,col="brown")
abline(h=20,lty=2,lwd=2,col="brown")
abline(coef=c(6,1),lty=2,lwd=2,col="red")
abline(coef=c(-6,1),lty=2,lwd=2,col="red")

par(mar=c(5,5,1,1))
plot(c$V1, c$V2, pch=20, cex=c$V4*5, xlab="Pavo muticus Indv 1", ylab="Pavo muticus Indv 3", font=2, font.lab=4, cex.lab=1.5, cex.main=2, cex.axis=1.5, col=ifelse(c$V4>quantile(c$V4,0.99), "red", ifelse(c$V4>quantile(c$V4,0.95), "blue", ifelse(c$V4>quantile(c$V4,0.9), "green", "black"))))
abline(v=60,lty=2,lwd=2,col="brown")
abline(h=60,lty=2,lwd=2,col="brown")
abline(v=40,lty=2,lwd=2,col="brown")
abline(h=40,lty=2,lwd=2,col="brown")
abline(v=20,lty=2,lwd=2,col="brown")
abline(h=20,lty=2,lwd=2,col="brown")
abline(coef=c(6,1),lty=2,lwd=2,col="red")
abline(coef=c(-6,1),lty=2,lwd=2,col="red")

par(mar=c(5,5,1,1))
plot(d$V1, d$V2, pch=20, cex=d$V4*5, xlab="Pavo muticus Indv 1", ylab="Pavo muticus Indv 4", font=2, font.lab=4, cex.lab=1.5, cex.main=2, cex.axis=1.5, col=ifelse(d$V4>quantile(d$V4,0.99), "red", ifelse(d$V4>quantile(d$V4,0.95), "blue", ifelse(d$V4>quantile(d$V4,0.9), "green", "black"))))
abline(v=60,lty=2,lwd=2,col="brown")
abline(h=60,lty=2,lwd=2,col="brown")
abline(v=40,lty=2,lwd=2,col="brown")
abline(h=40,lty=2,lwd=2,col="brown")
abline(v=20,lty=2,lwd=2,col="brown")
abline(h=20,lty=2,lwd=2,col="brown")
abline(coef=c(6,1),lty=2,lwd=2,col="red")
abline(coef=c(-6,1),lty=2,lwd=2,col="red")

dev.off()

###############################################################
#Figure S2
a <- read.table("pm1_pc_ai_distr.txt")
(a$V3*100)/sum(a$V3)->a$V4
b <- read.table("pm2_pc_ai_distr.txt")
(b$V3*100)/sum(b$V3)->b$V4
c <- read.table("pm3_pc_ai_distr.txt")
(c$V3*100)/sum(c$V3)->c$V4
d <- read.table("pm4_pc_ai_distr.txt")
(d$V3*100)/sum(d$V3)->d$V4


pdf("PC_4panel_distr_diff_ne.pdf",width=10,height=10)
par(mfrow=c(2,2))

par(mar=c(5,5,1,1))
plot(a$V1, a$V2, pch=20, cex=a$V4*5, xlab="Pavo muticus Indv 1", ylab="Pavo cristatus", font=2, font.lab=4, cex.lab=1.5, cex.main=2, cex.axis=1.5, col=ifelse(a$V4>quantile(a$V4,0.99), "red", ifelse(a$V4>quantile(a$V4,0.95), "blue", ifelse(a$V4>quantile(a$V4,0.9), "green", "black"))))
abline(v=60,lty=2,lwd=2,col="brown")
abline(h=60,lty=2,lwd=2,col="brown")
abline(v=40,lty=2,lwd=2,col="brown")
abline(h=40,lty=2,lwd=2,col="brown")
abline(v=20,lty=2,lwd=2,col="brown")
abline(h=20,lty=2,lwd=2,col="brown")
abline(coef=c(6,1),lty=2,lwd=2,col="red")
abline(coef=c(-6,1),lty=2,lwd=2,col="red")

par(mar=c(5,5,1,1))
plot(b$V1, b$V2, pch=20, cex=b$V4*5, xlab="Pavo muticus Indv 2", ylab="Pavo cristatus", font=2, font.lab=4, cex.lab=1.5, cex.main=2, cex.axis=1.5, col=ifelse(b$V4>quantile(b$V4,0.99), "red", ifelse(b$V4>quantile(b$V4,0.95), "blue", ifelse(b$V4>quantile(b$V4,0.9), "green", "black"))))
abline(v=60,lty=2,lwd=2,col="brown")
abline(h=60,lty=2,lwd=2,col="brown")
abline(v=40,lty=2,lwd=2,col="brown")
abline(h=40,lty=2,lwd=2,col="brown")
abline(v=20,lty=2,lwd=2,col="brown")
abline(h=20,lty=2,lwd=2,col="brown")
abline(coef=c(6,1),lty=2,lwd=2,col="red")
abline(coef=c(-6,1),lty=2,lwd=2,col="red")

par(mar=c(5,5,1,1))
plot(c$V1, c$V2, pch=20, cex=c$V4*5, xlab="Pavo muticus Indv 3", ylab="Pavo cristatus", font=2, font.lab=4, cex.lab=1.5, cex.main=2, cex.axis=1.5, col=ifelse(c$V4>quantile(c$V4,0.99), "red", ifelse(c$V4>quantile(c$V4,0.95), "blue", ifelse(c$V4>quantile(c$V4,0.9), "green", "black"))))
abline(v=60,lty=2,lwd=2,col="brown")
abline(h=60,lty=2,lwd=2,col="brown")
abline(v=40,lty=2,lwd=2,col="brown")
abline(h=40,lty=2,lwd=2,col="brown")
abline(v=20,lty=2,lwd=2,col="brown")
abline(h=20,lty=2,lwd=2,col="brown")
abline(coef=c(6,1),lty=2,lwd=2,col="red")
abline(coef=c(-6,1),lty=2,lwd=2,col="red")

par(mar=c(5,5,1,1))
plot(d$V1, d$V2, pch=20, cex=d$V4*5, xlab="Pavo muticus Indv 4", ylab="Pavo cristatus", font=2, font.lab=4, cex.lab=1.5, cex.main=2, cex.axis=1.5, col=ifelse(d$V4>quantile(d$V4,0.99), "red", ifelse(d$V4>quantile(d$V4,0.95), "blue", ifelse(d$V4>quantile(d$V4,0.9), "green", "black"))))
abline(v=60,lty=2,lwd=2,col="brown")
abline(h=60,lty=2,lwd=2,col="brown")
abline(v=40,lty=2,lwd=2,col="brown")
abline(h=40,lty=2,lwd=2,col="brown")
abline(v=20,lty=2,lwd=2,col="brown")
abline(h=20,lty=2,lwd=2,col="brown")
abline(coef=c(6,1),lty=2,lwd=2,col="red")
abline(coef=c(-6,1),lty=2,lwd=2,col="red")

dev.off()
