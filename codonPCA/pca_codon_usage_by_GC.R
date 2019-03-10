# Codon usage bias across kingdoms divided by its GC content 
# Collaboration with Olivier Jaillon for TARA project

setwd("/Users/enovoa/kellislab/PROJECTS/CODON_USAGE_ACROSS_SPECIES_GC_CONTENT/RESULTS/EMBL_CDS_ALL_SEQUENCES_1625_SPS/results_by_species/all_joint")

# 1. GET DATA

data_allGC<-read.table("all_allGC.txt.transposed", header=T)
data_below30<-read.table("all_below_30.txt.transposed", header=T)
data_30_40<-read.table("all_30-40.txt.transposed", header=T)
data_40_50<-read.table("all_40-50.txt.transposed", header=T)
data_50_60<-read.table("all_50-60.txt.transposed", header=T)
data_60_70<-read.table("all_60-70.txt.transposed", header=T)

dat_allGC<-data_allGC[,2:64]
dat_below30<-data_below30[,2:64]
dat_30_40<-data_30_40[,2:64]
dat_40_50<-data_40_50[,2:64]
dat_50_60<-data_50_60[,2:64]
dat_60_70<-data_60_70[,2:64]

#rownames(dat_allGC)<-data_allGC$Species
rownames(dat_below30)<-data_below30$Species
rownames(dat_30_40)<-data_30_40$Species
rownames(dat_40_50)<-data_40_50$Species
rownames(dat_50_60)<-data_50_60$Species
rownames(dat_60_70)<-data_60_70$Species


# 2. BUILD PCA
model_allGC<-prcomp(dat_allGC)
model_below30<-prcomp(dat_below30)
model_30_40<-prcomp(dat_30_40)
model_40_50<-prcomp(dat_40_50)
model_50_60<-prcomp(dat_50_60)
model_60_70<-prcomp(dat_60_70)

#summary(model_allGC)
#summary(model_below30)
#summary(model_30_40)
#summary(model_40_50)
#summary(model_50_60)
#summary(model_60_70)


# 3. GET SCORES MATRIX PCA 
scores_allGC<-model_allGC$x
scores_below30<-model_below30$x
scores_30_40<-model_30_40$x
scores_40_50<-model_40_50$x
scores_50_60<-model_50_60$x
scores_60_70<-model_60_70$x

scores_data_allGC<-as.matrix(scores_allGC[,1:3])
scores_data_below30<-as.matrix(scores_below30[,1:3])
scores_data_30_40<-as.matrix(scores_30_40[,1:3])
scores_data_40_50<-as.matrix(scores_40_50[,1:3])
scores_data_50_60<-as.matrix(scores_50_60[,1:3])
scores_data_60_70<-as.matrix(scores_60_70[,1:3])

#rownames(scores_data_allGC)
#rownames(scores_data_below30)
#rownames(scores_data_30_40)
#rownames(scores_data_40_50)
#rownames(scores_data_50_60)
#rownames(scores_data_60_70)

pc1_all_allGC<-as.vector(scores_data_allGC[,1])
pc2_all_allGC<-as.vector(scores_data_allGC[,2])
pc3_all_allGC<-as.vector(scores_data_allGC[,3])
pc1_all_below30<-as.vector(scores_data_below30[,1])
pc2_all_below30<-as.vector(scores_data_below30[,2])
pc3_all_below30<-as.vector(scores_data_below30[,3])
pc1_all_30_40<-as.vector(scores_data_30_40[,1])
pc2_all_30_40<-as.vector(scores_data_30_40[,2])
pc3_all_30_40<-as.vector(scores_data_30_40[,3])
pc1_all_40_50<-as.vector(scores_data_40_50[,1])
pc2_all_40_50<-as.vector(scores_data_40_50[,2])
pc3_all_40_50<-as.vector(scores_data_40_50[,3])
pc1_all_50_60<-as.vector(scores_data_50_60[,1])
pc2_all_50_60<-as.vector(scores_data_50_60[,2])
pc3_all_50_60<-as.vector(scores_data_50_60[,3])
pc1_all_60_70<-as.vector(scores_data_60_70[,1])
pc2_all_60_70<-as.vector(scores_data_60_70[,2])
pc3_all_60_70<-as.vector(scores_data_60_70[,3])


# 4. GET SCORES FOR EUK AND BACT

arch_allGC<-scores_data_allGC[1:120,] 		# 120
bact_allGC<-scores_data_allGC[121:1461,]  	# 1341
euk_allGC<-scores_data_allGC[1462:1625,]	# 164 --> total = 1625
pc1_arch_allGC<-as.vector(arch_allGC[,1])
pc2_arch_allGC<-as.vector(arch_allGC[,2])
pc3_arch_allGC<-as.vector(arch_allGC[,3])
pc1_euk_allGC<-as.vector(euk_allGC[,1])
pc2_euk_allGC<-as.vector(euk_allGC[,2])
pc3_euk_allGC<-as.vector(euk_allGC[,3])
pc1_bact_allGC<-as.vector(bact_allGC[,1])
pc2_bact_allGC<-as.vector(bact_allGC[,2])
pc3_bact_allGC<-as.vector(bact_allGC[,3])

#arch_below30<-scores_data_below30[1:1,] 	# 1
bact_below30<-scores_data_below30[2:37,]  	# 36
euk_below30<-scores_data_below30[38:47,]	# 10 --> total = 47
#pc1_arch_below30<-as.vector(arch_below30[,1])
#pc2_arch_below30<-as.vector(arch_below30[,2])
#pc3_arch_below30<-as.vector(arch_below30[,3])
pc1_euk_below30<-as.vector(euk_below30[,1])
pc2_euk_below30<-as.vector(euk_below30[,2])
pc3_euk_below30<-as.vector(euk_below30[,3])
pc1_bact_below30<-as.vector(bact_below30[,1])
pc2_bact_below30<-as.vector(bact_below30[,2])
pc3_bact_below30<-as.vector(bact_below30[,3])

arch_30_40<-scores_data_30_40[1:26,]	# 26
bact_30_40<-scores_data_30_40[27:331,]  # 305
euk_30_40<-scores_data_30_40[332:356,]	# 25 --> total = 356
pc1_arch_30_40<-as.vector(arch_30_40[,1])
pc2_arch_30_40<-as.vector(arch_30_40[,2])
pc3_arch_30_40<-as.vector(arch_30_40[,3])
pc1_euk_30_40<-as.vector(euk_30_40[,1])
pc2_euk_30_40<-as.vector(euk_30_40[,2])
pc3_euk_30_40<-as.vector(euk_30_40[,3])
pc1_bact_30_40<-as.vector(bact_30_40[,1])
pc2_bact_30_40<-as.vector(bact_30_40[,2])
pc3_bact_30_40<-as.vector(bact_30_40[,3])

arch_40_50<-scores_data_40_50[1:35,]	# 35
bact_40_50<-scores_data_40_50[36:375,]  # 340
euk_40_50<-scores_data_40_50[376:439,]	# 64 --> total = 439
pc1_arch_40_50<-as.vector(arch_40_50[,1])
pc2_arch_40_50<-as.vector(arch_40_50[,2])
pc3_arch_40_50<-as.vector(arch_40_50[,3])
pc1_euk_40_50<-as.vector(euk_40_50[,1])
pc2_euk_40_50<-as.vector(euk_40_50[,2])
pc3_euk_40_50<-as.vector(euk_40_50[,3])
pc1_bact_40_50<-as.vector(bact_40_50[,1])
pc2_bact_40_50<-as.vector(bact_40_50[,2])
pc3_bact_40_50<-as.vector(bact_40_50[,3])

arch_50_60<-scores_data_50_60[1:34,]	# 34
bact_50_60<-scores_data_50_60[35:327,]  # 293
euk_50_60<-scores_data_50_60[328:377,]	# 50 --> total = 377
pc1_arch_50_60<-as.vector(arch_50_60[,1])
pc2_arch_50_60<-as.vector(arch_50_60[,2])
pc3_arch_50_60<-as.vector(arch_50_60[,3])
pc1_euk_50_60<-as.vector(euk_50_60[,1])
pc2_euk_50_60<-as.vector(euk_50_60[,2])
pc3_euk_50_60<-as.vector(euk_50_60[,3])
pc1_bact_50_60<-as.vector(bact_50_60[,1])
pc2_bact_50_60<-as.vector(bact_50_60[,2])
pc3_bact_50_60<-as.vector(bact_50_60[,3])

arch_60_70<-scores_data_60_70[1:24,]	# 24
bact_60_70<-scores_data_60_70[25:326,]  # 302
euk_60_70<-scores_data_60_70[327:340,]	# 14 --> total = 340
pc1_arch_60_70<-as.vector(arch_60_70[,1])
pc2_arch_60_70<-as.vector(arch_60_70[,2])
pc3_arch_60_70<-as.vector(arch_60_70[,3])
pc1_euk_60_70<-as.vector(euk_60_70[,1])
pc2_euk_60_70<-as.vector(euk_60_70[,2])
pc3_euk_60_70<-as.vector(euk_60_70[,3])
pc1_bact_60_70<-as.vector(bact_60_70[,1])
pc2_bact_60_70<-as.vector(bact_60_70[,2])
pc3_bact_60_70<-as.vector(bact_60_70[,3])

# 5. GET LOADINGS MATRIX PCA 
loadings_allGC<-model_allGC$rotation
loadings_below30<-model_below30$rotation
loadings_30_40<-model_30_40$rotation
loadings_40_50<-model_40_50$rotation
loadings_50_60<-model_50_60$rotation
loadings_60_70<-model_60_70$rotation

loadings_data_allGC<-as.matrix(loadings_allGC[,1:3])
loadings_data_below30<-as.matrix(loadings_below30[,1:3])
loadings_data_30_40<-as.matrix(loadings_30_40[,1:3])
loadings_data_40_50<-as.matrix(loadings_40_50[,1:3])
loadings_data_50_60<-as.matrix(loadings_50_60[,1:3])
loadings_data_60_70<-as.matrix(loadings_60_70[,1:3])

loadings_pc1_all_allGC<-as.vector(loadings_data_allGC[,1])
loadings_pc2_all_allGC<-as.vector(loadings_data_allGC[,2])
loadings_pc3_all_allGC<-as.vector(loadings_data_allGC[,3])
loadings_pc1_all_below30<-as.vector(loadings_data_below30[,1])
loadings_pc2_all_below30<-as.vector(loadings_data_below30[,2])
loadings_pc3_all_below30<-as.vector(loadings_data_below30[,3])
loadings_pc1_all_30_40<-as.vector(loadings_data_30_40[,1])
loadings_pc2_all_30_40<-as.vector(loadings_data_30_40[,2])
loadings_pc3_all_30_40<-as.vector(loadings_data_30_40[,3])
loadings_pc1_all_40_50<-as.vector(loadings_data_40_50[,1])
loadings_pc2_all_40_50<-as.vector(loadings_data_40_50[,2])
loadings_pc3_all_40_50<-as.vector(loadings_data_40_50[,3])
loadings_pc1_all_50_60<-as.vector(loadings_data_50_60[,1])
loadings_pc2_all_50_60<-as.vector(loadings_data_50_60[,2])
loadings_pc3_all_50_60<-as.vector(loadings_data_50_60[,3])
loadings_pc1_all_60_70<-as.vector(loadings_data_60_70[,1])
loadings_pc2_all_60_70<-as.vector(loadings_data_60_70[,2])
loadings_pc3_all_60_70<-as.vector(loadings_data_60_70[,3])


# PLOT PCA SCORES PC1-PC2

pdf(file="scores_by_kingdom_RGF_allGC.pdf", height=7, width=7)
plot(pc1_all_allGC, pc2_all_allGC, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (all GC's)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc1_arch_allGC, pc2_arch_allGC, pch=15, col="red")
points(pc1_bact_allGC, pc2_bact_allGC, pch=16, col="purple")
points(pc1_euk_allGC, pc2_euk_allGC, pch=17,col="forestgreen")
#arrows(0, 0, X[,1], X[,2], len=0.1, col="black")
#text(1.4*X, rownames(X), col="black", xpd=T)
legend('topright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="scores_by_kingdom_RGF_below30.pdf", height=7, width=7)
plot(pc1_all_below30, pc2_all_below30, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (GC < 30%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
#points(pc1_arch_below30, pc2_arch_below30, pch=15, col="red")
points(pc1_bact_below30, pc2_bact_below30, pch=16, col="purple")
points(pc1_euk_below30, pc2_euk_below30, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Bacteria','Eukarya'),col=c('purple','forestgreen'),pch=c(16,17), box.col="white")
dev.off()

pdf(file="scores_by_kingdom_RGF_30_40.pdf", height=7, width=7)
plot(pc1_all_30_40, pc2_all_30_40, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (GC = 30-40%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc1_arch_30_40, pc2_arch_30_40, pch=15, col="red")
points(pc1_bact_30_40, pc2_bact_30_40, pch=16, col="purple")
points(pc1_euk_30_40, pc2_euk_30_40, pch=17,col="forestgreen")
legend('bottomright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="scores_by_kingdom_RGF_40_50.pdf", height=7, width=7)
plot(pc1_all_40_50, pc2_all_40_50, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (GC = 40-50%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc1_arch_40_50, pc2_arch_40_50, pch=15, col="red")
points(pc1_bact_40_50, pc2_bact_40_50, pch=16, col="purple")
points(pc1_euk_40_50, pc2_euk_40_50, pch=17,col="forestgreen")
legend('bottomright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="scores_by_kingdom_RGF_50_60.pdf", height=7, width=7)
plot(pc1_all_50_60, pc2_all_50_60, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (GC = 50-60%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc1_arch_50_60, pc2_arch_50_60, pch=15, col="red")
points(pc1_bact_50_60, pc2_bact_50_60, pch=16, col="purple")
points(pc1_euk_50_60, pc2_euk_50_60, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="scores_by_kingdom_RGF_60_70.pdf", height=7, width=7)
plot(pc1_all_60_70, pc2_all_60_70, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (GC = 60-70%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc1_arch_60_70, pc2_arch_60_70, pch=15, col="red")
points(pc1_bact_60_70, pc2_bact_60_70, pch=16, col="purple")
points(pc1_euk_60_70, pc2_euk_60_70, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()


# PLOT PCA LOADINGS PC1-PC2

# Separate labels in order to color according to its last codon --> use data_allGC, which can be used for all plots
labels_A<-cbind(rownames(loadings_data_allGC)[1],rownames(loadings_data_allGC)[5],rownames(loadings_data_allGC)[9],rownames(loadings_data_allGC)[13],rownames(loadings_data_allGC)[17],rownames(loadings_data_allGC)[21],rownames(loadings_data_allGC)[25],rownames(loadings_data_allGC)[29],rownames(loadings_data_allGC)[33],rownames(loadings_data_allGC)[37],rownames(loadings_data_allGC)[41],rownames(loadings_data_allGC)[45],rownames(loadings_data_allGC)[49],rownames(loadings_data_allGC)[53],rownames(loadings_data_allGC)[57],rownames(loadings_data_allGC)[61])
labels_C<-cbind(rownames(loadings_data_allGC)[2],rownames(loadings_data_allGC)[6],rownames(loadings_data_allGC)[10],rownames(loadings_data_allGC)[14],rownames(loadings_data_allGC)[18],rownames(loadings_data_allGC)[22],rownames(loadings_data_allGC)[26],rownames(loadings_data_allGC)[30],rownames(loadings_data_allGC)[34],rownames(loadings_data_allGC)[38],rownames(loadings_data_allGC)[42],rownames(loadings_data_allGC)[46],rownames(loadings_data_allGC)[50],rownames(loadings_data_allGC)[54],rownames(loadings_data_allGC)[58],rownames(loadings_data_allGC)[62])
labels_G<-cbind(rownames(loadings_data_allGC)[3],rownames(loadings_data_allGC)[7],rownames(loadings_data_allGC)[11],rownames(loadings_data_allGC)[15],rownames(loadings_data_allGC)[19],rownames(loadings_data_allGC)[23],rownames(loadings_data_allGC)[27],rownames(loadings_data_allGC)[31],rownames(loadings_data_allGC)[35],rownames(loadings_data_allGC)[39],rownames(loadings_data_allGC)[43],rownames(loadings_data_allGC)[47],rownames(loadings_data_allGC)[51],rownames(loadings_data_allGC)[55],rownames(loadings_data_allGC)[59],rownames(loadings_data_allGC)[63])
labels_U<-cbind(rownames(loadings_data_allGC)[4],rownames(loadings_data_allGC)[8],rownames(loadings_data_allGC)[12],rownames(loadings_data_allGC)[16],rownames(loadings_data_allGC)[20],rownames(loadings_data_allGC)[24],rownames(loadings_data_allGC)[28],rownames(loadings_data_allGC)[32],rownames(loadings_data_allGC)[36],rownames(loadings_data_allGC)[40],rownames(loadings_data_allGC)[44],rownames(loadings_data_allGC)[48],rownames(loadings_data_allGC)[52],rownames(loadings_data_allGC)[56],rownames(loadings_data_allGC)[60],rownames(loadings_data_allGC)[64])



# a) all GC
pdf(file="loadings_by_kingdom_RGF_allGC.pdf", height=7, width=7)
plot(loadings_pc1_all_allGC,loadings_pc2_all_allGC, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (all GC's)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)

# loadings pc1
loadings_pc1_all_allGC_A<-cbind(loadings_pc1_all_allGC[1],loadings_pc1_all_allGC[5],loadings_pc1_all_allGC[9],loadings_pc1_all_allGC[13],loadings_pc1_all_allGC[17],loadings_pc1_all_allGC[21],loadings_pc1_all_allGC[25],loadings_pc1_all_allGC[29],loadings_pc1_all_allGC[33],loadings_pc1_all_allGC[37],loadings_pc1_all_allGC[41],loadings_pc1_all_allGC[45],loadings_pc1_all_allGC[49],loadings_pc1_all_allGC[53],loadings_pc1_all_allGC[57],loadings_pc1_all_allGC[61])

loadings_pc1_all_allGC_C<-cbind(loadings_pc1_all_allGC[2],loadings_pc1_all_allGC[6],loadings_pc1_all_allGC[10],loadings_pc1_all_allGC[14],loadings_pc1_all_allGC[18],loadings_pc1_all_allGC[22],loadings_pc1_all_allGC[26],loadings_pc1_all_allGC[30],loadings_pc1_all_allGC[34],loadings_pc1_all_allGC[38],loadings_pc1_all_allGC[42],loadings_pc1_all_allGC[46],loadings_pc1_all_allGC[50],loadings_pc1_all_allGC[54],loadings_pc1_all_allGC[58],loadings_pc1_all_allGC[62])

loadings_pc1_all_allGC_G<-cbind(loadings_pc1_all_allGC[3],loadings_pc1_all_allGC[7],loadings_pc1_all_allGC[11],loadings_pc1_all_allGC[15],loadings_pc1_all_allGC[19],loadings_pc1_all_allGC[23],loadings_pc1_all_allGC[27],loadings_pc1_all_allGC[31],loadings_pc1_all_allGC[35],loadings_pc1_all_allGC[39],loadings_pc1_all_allGC[43],loadings_pc1_all_allGC[47],loadings_pc1_all_allGC[51],loadings_pc1_all_allGC[55],loadings_pc1_all_allGC[59],loadings_pc1_all_allGC[63])

loadings_pc1_all_allGC_U<-cbind(loadings_pc1_all_allGC[4],loadings_pc1_all_allGC[8],loadings_pc1_all_allGC[12],loadings_pc1_all_allGC[16],loadings_pc1_all_allGC[20],loadings_pc1_all_allGC[24],loadings_pc1_all_allGC[28],loadings_pc1_all_allGC[32],loadings_pc1_all_allGC[36],loadings_pc1_all_allGC[40],loadings_pc1_all_allGC[44],loadings_pc1_all_allGC[48],loadings_pc1_all_allGC[52],loadings_pc1_all_allGC[56],loadings_pc1_all_allGC[60],loadings_pc1_all_allGC[64])

# loadings pc2

loadings_pc2_all_allGC_A<-cbind(loadings_pc2_all_allGC[1],loadings_pc2_all_allGC[5],loadings_pc2_all_allGC[9],loadings_pc2_all_allGC[13],loadings_pc2_all_allGC[17],loadings_pc2_all_allGC[21],loadings_pc2_all_allGC[25],loadings_pc2_all_allGC[29],loadings_pc2_all_allGC[33],loadings_pc2_all_allGC[37],loadings_pc2_all_allGC[41],loadings_pc2_all_allGC[45],loadings_pc2_all_allGC[49],loadings_pc2_all_allGC[53],loadings_pc2_all_allGC[57],loadings_pc2_all_allGC[61])

loadings_pc2_all_allGC_C<-cbind(loadings_pc2_all_allGC[2],loadings_pc2_all_allGC[6],loadings_pc2_all_allGC[10],loadings_pc2_all_allGC[14],loadings_pc2_all_allGC[18],loadings_pc2_all_allGC[22],loadings_pc2_all_allGC[26],loadings_pc2_all_allGC[30],loadings_pc2_all_allGC[34],loadings_pc2_all_allGC[38],loadings_pc2_all_allGC[42],loadings_pc2_all_allGC[46],loadings_pc2_all_allGC[50],loadings_pc2_all_allGC[54],loadings_pc2_all_allGC[58],loadings_pc2_all_allGC[62])

loadings_pc2_all_allGC_G<-cbind(loadings_pc2_all_allGC[3],loadings_pc2_all_allGC[7],loadings_pc2_all_allGC[11],loadings_pc2_all_allGC[15],loadings_pc2_all_allGC[19],loadings_pc2_all_allGC[23],loadings_pc2_all_allGC[27],loadings_pc2_all_allGC[31],loadings_pc2_all_allGC[35],loadings_pc2_all_allGC[39],loadings_pc2_all_allGC[43],loadings_pc2_all_allGC[47],loadings_pc2_all_allGC[51],loadings_pc2_all_allGC[55],loadings_pc2_all_allGC[59],loadings_pc2_all_allGC[63])

loadings_pc2_all_allGC_U<-cbind(loadings_pc2_all_allGC[4],loadings_pc2_all_allGC[8],loadings_pc2_all_allGC[12],loadings_pc2_all_allGC[16],loadings_pc2_all_allGC[20],loadings_pc2_all_allGC[24],loadings_pc2_all_allGC[28],loadings_pc2_all_allGC[32],loadings_pc2_all_allGC[36],loadings_pc2_all_allGC[40],loadings_pc2_all_allGC[44],loadings_pc2_all_allGC[48],loadings_pc2_all_allGC[52],loadings_pc2_all_allGC[56],loadings_pc2_all_allGC[60],loadings_pc2_all_allGC[64])


text(loadings_pc1_all_allGC_A,loadings_pc2_all_allGC_A, labels_A, col="blue", xpd=T)
text(loadings_pc1_all_allGC_C,loadings_pc2_all_allGC_C, labels_C, col="red", xpd=T)
text(loadings_pc1_all_allGC_G,loadings_pc2_all_allGC_G, labels_G, col="orange", xpd=T)
text(loadings_pc1_all_allGC_U,loadings_pc2_all_allGC_U, labels_U, col="purple", xpd=T)
legend('bottomright', cex=1,ncol=1,c('A-ended','U-ended','G-ended','C-ended'),col=c('blue','purple','orange','red'),pch=c(15,15,15,15), box.col="white")
dev.off()

# b) Below 30
pdf(file="loadings_by_kingdom_RGF_below30.pdf", height=7, width=7)
plot(loadings_pc1_all_below30,loadings_pc2_all_below30, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (GC < 30%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)

# loadings pc1
loadings_pc1_all_below30_A<-cbind(loadings_pc1_all_below30[1],loadings_pc1_all_below30[5],loadings_pc1_all_below30[9],loadings_pc1_all_below30[13],loadings_pc1_all_below30[17],loadings_pc1_all_below30[21],loadings_pc1_all_below30[25],loadings_pc1_all_below30[29],loadings_pc1_all_below30[33],loadings_pc1_all_below30[37],loadings_pc1_all_below30[41],loadings_pc1_all_below30[45],loadings_pc1_all_below30[49],loadings_pc1_all_below30[53],loadings_pc1_all_below30[57],loadings_pc1_all_below30[61])

loadings_pc1_all_below30_C<-cbind(loadings_pc1_all_below30[2],loadings_pc1_all_below30[6],loadings_pc1_all_below30[10],loadings_pc1_all_below30[14],loadings_pc1_all_below30[18],loadings_pc1_all_below30[22],loadings_pc1_all_below30[26],loadings_pc1_all_below30[30],loadings_pc1_all_below30[34],loadings_pc1_all_below30[38],loadings_pc1_all_below30[42],loadings_pc1_all_below30[46],loadings_pc1_all_below30[50],loadings_pc1_all_below30[54],loadings_pc1_all_below30[58],loadings_pc1_all_below30[62])

loadings_pc1_all_below30_G<-cbind(loadings_pc1_all_below30[3],loadings_pc1_all_below30[7],loadings_pc1_all_below30[11],loadings_pc1_all_below30[15],loadings_pc1_all_below30[19],loadings_pc1_all_below30[23],loadings_pc1_all_below30[27],loadings_pc1_all_below30[31],loadings_pc1_all_below30[35],loadings_pc1_all_below30[39],loadings_pc1_all_below30[43],loadings_pc1_all_below30[47],loadings_pc1_all_below30[51],loadings_pc1_all_below30[55],loadings_pc1_all_below30[59],loadings_pc1_all_below30[63])

loadings_pc1_all_below30_U<-cbind(loadings_pc1_all_below30[4],loadings_pc1_all_below30[8],loadings_pc1_all_below30[12],loadings_pc1_all_below30[16],loadings_pc1_all_below30[20],loadings_pc1_all_below30[24],loadings_pc1_all_below30[28],loadings_pc1_all_below30[32],loadings_pc1_all_below30[36],loadings_pc1_all_below30[40],loadings_pc1_all_below30[44],loadings_pc1_all_below30[48],loadings_pc1_all_below30[52],loadings_pc1_all_below30[56],loadings_pc1_all_below30[60],loadings_pc1_all_below30[64])

# loadings pc2

loadings_pc2_all_below30_A<-cbind(loadings_pc2_all_below30[1],loadings_pc2_all_below30[5],loadings_pc2_all_below30[9],loadings_pc2_all_below30[13],loadings_pc2_all_below30[17],loadings_pc2_all_below30[21],loadings_pc2_all_below30[25],loadings_pc2_all_below30[29],loadings_pc2_all_below30[33],loadings_pc2_all_below30[37],loadings_pc2_all_below30[41],loadings_pc2_all_below30[45],loadings_pc2_all_below30[49],loadings_pc2_all_below30[53],loadings_pc2_all_below30[57],loadings_pc2_all_below30[61])

loadings_pc2_all_below30_C<-cbind(loadings_pc2_all_below30[2],loadings_pc2_all_below30[6],loadings_pc2_all_below30[10],loadings_pc2_all_below30[14],loadings_pc2_all_below30[18],loadings_pc2_all_below30[22],loadings_pc2_all_below30[26],loadings_pc2_all_below30[30],loadings_pc2_all_below30[34],loadings_pc2_all_below30[38],loadings_pc2_all_below30[42],loadings_pc2_all_below30[46],loadings_pc2_all_below30[50],loadings_pc2_all_below30[54],loadings_pc2_all_below30[58],loadings_pc2_all_below30[62])

loadings_pc2_all_below30_G<-cbind(loadings_pc2_all_below30[3],loadings_pc2_all_below30[7],loadings_pc2_all_below30[11],loadings_pc2_all_below30[15],loadings_pc2_all_below30[19],loadings_pc2_all_below30[23],loadings_pc2_all_below30[27],loadings_pc2_all_below30[31],loadings_pc2_all_below30[35],loadings_pc2_all_below30[39],loadings_pc2_all_below30[43],loadings_pc2_all_below30[47],loadings_pc2_all_below30[51],loadings_pc2_all_below30[55],loadings_pc2_all_below30[59],loadings_pc2_all_below30[63])

loadings_pc2_all_below30_U<-cbind(loadings_pc2_all_below30[4],loadings_pc2_all_below30[8],loadings_pc2_all_below30[12],loadings_pc2_all_below30[16],loadings_pc2_all_below30[20],loadings_pc2_all_below30[24],loadings_pc2_all_below30[28],loadings_pc2_all_below30[32],loadings_pc2_all_below30[36],loadings_pc2_all_below30[40],loadings_pc2_all_below30[44],loadings_pc2_all_below30[48],loadings_pc2_all_below30[52],loadings_pc2_all_below30[56],loadings_pc2_all_below30[60],loadings_pc2_all_below30[64])


text(loadings_pc1_all_below30_A,loadings_pc2_all_below30_A, labels_A, col="blue", xpd=T)
text(loadings_pc1_all_below30_C,loadings_pc2_all_below30_C, labels_C, col="red", xpd=T)
text(loadings_pc1_all_below30_G,loadings_pc2_all_below30_G, labels_G, col="orange", xpd=T)
text(loadings_pc1_all_below30_U,loadings_pc2_all_below30_U, labels_U, col="purple", xpd=T)
legend('bottomright', cex=1,ncol=1,c('A-ended','U-ended','G-ended','C-ended'),col=c('blue','purple','orange','red'),pch=c(15,15,15,15), box.col="white")
dev.off()


# c) 30-40%
pdf(file="loadings_by_kingdom_RGF_30_40.pdf", height=7, width=7)
plot(loadings_pc1_all_30_40,loadings_pc2_all_30_40, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (GC = 30-40%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)

# loadings pc1
loadings_pc1_all_30_40_A<-cbind(loadings_pc1_all_30_40[1],loadings_pc1_all_30_40[5],loadings_pc1_all_30_40[9],loadings_pc1_all_30_40[13],loadings_pc1_all_30_40[17],loadings_pc1_all_30_40[21],loadings_pc1_all_30_40[25],loadings_pc1_all_30_40[29],loadings_pc1_all_30_40[33],loadings_pc1_all_30_40[37],loadings_pc1_all_30_40[41],loadings_pc1_all_30_40[45],loadings_pc1_all_30_40[49],loadings_pc1_all_30_40[53],loadings_pc1_all_30_40[57],loadings_pc1_all_30_40[61])

loadings_pc1_all_30_40_C<-cbind(loadings_pc1_all_30_40[2],loadings_pc1_all_30_40[6],loadings_pc1_all_30_40[10],loadings_pc1_all_30_40[14],loadings_pc1_all_30_40[18],loadings_pc1_all_30_40[22],loadings_pc1_all_30_40[26],loadings_pc1_all_30_40[30],loadings_pc1_all_30_40[34],loadings_pc1_all_30_40[38],loadings_pc1_all_30_40[42],loadings_pc1_all_30_40[46],loadings_pc1_all_30_40[50],loadings_pc1_all_30_40[54],loadings_pc1_all_30_40[58],loadings_pc1_all_30_40[62])

loadings_pc1_all_30_40_G<-cbind(loadings_pc1_all_30_40[3],loadings_pc1_all_30_40[7],loadings_pc1_all_30_40[11],loadings_pc1_all_30_40[15],loadings_pc1_all_30_40[19],loadings_pc1_all_30_40[23],loadings_pc1_all_30_40[27],loadings_pc1_all_30_40[31],loadings_pc1_all_30_40[35],loadings_pc1_all_30_40[39],loadings_pc1_all_30_40[43],loadings_pc1_all_30_40[47],loadings_pc1_all_30_40[51],loadings_pc1_all_30_40[55],loadings_pc1_all_30_40[59],loadings_pc1_all_30_40[63])

loadings_pc1_all_30_40_U<-cbind(loadings_pc1_all_30_40[4],loadings_pc1_all_30_40[8],loadings_pc1_all_30_40[12],loadings_pc1_all_30_40[16],loadings_pc1_all_30_40[20],loadings_pc1_all_30_40[24],loadings_pc1_all_30_40[28],loadings_pc1_all_30_40[32],loadings_pc1_all_30_40[36],loadings_pc1_all_30_40[40],loadings_pc1_all_30_40[44],loadings_pc1_all_30_40[48],loadings_pc1_all_30_40[52],loadings_pc1_all_30_40[56],loadings_pc1_all_30_40[60],loadings_pc1_all_30_40[64])

# loadings pc2

loadings_pc2_all_30_40_A<-cbind(loadings_pc2_all_30_40[1],loadings_pc2_all_30_40[5],loadings_pc2_all_30_40[9],loadings_pc2_all_30_40[13],loadings_pc2_all_30_40[17],loadings_pc2_all_30_40[21],loadings_pc2_all_30_40[25],loadings_pc2_all_30_40[29],loadings_pc2_all_30_40[33],loadings_pc2_all_30_40[37],loadings_pc2_all_30_40[41],loadings_pc2_all_30_40[45],loadings_pc2_all_30_40[49],loadings_pc2_all_30_40[53],loadings_pc2_all_30_40[57],loadings_pc2_all_30_40[61])

loadings_pc2_all_30_40_C<-cbind(loadings_pc2_all_30_40[2],loadings_pc2_all_30_40[6],loadings_pc2_all_30_40[10],loadings_pc2_all_30_40[14],loadings_pc2_all_30_40[18],loadings_pc2_all_30_40[22],loadings_pc2_all_30_40[26],loadings_pc2_all_30_40[30],loadings_pc2_all_30_40[34],loadings_pc2_all_30_40[38],loadings_pc2_all_30_40[42],loadings_pc2_all_30_40[46],loadings_pc2_all_30_40[50],loadings_pc2_all_30_40[54],loadings_pc2_all_30_40[58],loadings_pc2_all_30_40[62])

loadings_pc2_all_30_40_G<-cbind(loadings_pc2_all_30_40[3],loadings_pc2_all_30_40[7],loadings_pc2_all_30_40[11],loadings_pc2_all_30_40[15],loadings_pc2_all_30_40[19],loadings_pc2_all_30_40[23],loadings_pc2_all_30_40[27],loadings_pc2_all_30_40[31],loadings_pc2_all_30_40[35],loadings_pc2_all_30_40[39],loadings_pc2_all_30_40[43],loadings_pc2_all_30_40[47],loadings_pc2_all_30_40[51],loadings_pc2_all_30_40[55],loadings_pc2_all_30_40[59],loadings_pc2_all_30_40[63])

loadings_pc2_all_30_40_U<-cbind(loadings_pc2_all_30_40[4],loadings_pc2_all_30_40[8],loadings_pc2_all_30_40[12],loadings_pc2_all_30_40[16],loadings_pc2_all_30_40[20],loadings_pc2_all_30_40[24],loadings_pc2_all_30_40[28],loadings_pc2_all_30_40[32],loadings_pc2_all_30_40[36],loadings_pc2_all_30_40[40],loadings_pc2_all_30_40[44],loadings_pc2_all_30_40[48],loadings_pc2_all_30_40[52],loadings_pc2_all_30_40[56],loadings_pc2_all_30_40[60],loadings_pc2_all_30_40[64])


text(loadings_pc1_all_30_40_A,loadings_pc2_all_30_40_A, labels_A, col="blue", xpd=T)
text(loadings_pc1_all_30_40_C,loadings_pc2_all_30_40_C, labels_C, col="red", xpd=T)
text(loadings_pc1_all_30_40_G,loadings_pc2_all_30_40_G, labels_G, col="orange", xpd=T)
text(loadings_pc1_all_30_40_U,loadings_pc2_all_30_40_U, labels_U, col="purple", xpd=T)
legend('bottomright', cex=1,ncol=1,c('A-ended','U-ended','G-ended','C-ended'),col=c('blue','purple','orange','red'),pch=c(15,15,15,15), box.col="white")
dev.off()

# d) 40-50%
pdf(file="loadings_by_kingdom_RGF_40_50.pdf", height=7, width=7)
plot(loadings_pc1_all_40_50,loadings_pc2_all_40_50, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (GC = 40-50%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)

# loadings pc1
loadings_pc1_all_40_50_A<-cbind(loadings_pc1_all_40_50[1],loadings_pc1_all_40_50[5],loadings_pc1_all_40_50[9],loadings_pc1_all_40_50[13],loadings_pc1_all_40_50[17],loadings_pc1_all_40_50[21],loadings_pc1_all_40_50[25],loadings_pc1_all_40_50[29],loadings_pc1_all_40_50[33],loadings_pc1_all_40_50[37],loadings_pc1_all_40_50[41],loadings_pc1_all_40_50[45],loadings_pc1_all_40_50[49],loadings_pc1_all_40_50[53],loadings_pc1_all_40_50[57],loadings_pc1_all_40_50[61])

loadings_pc1_all_40_50_C<-cbind(loadings_pc1_all_40_50[2],loadings_pc1_all_40_50[6],loadings_pc1_all_40_50[10],loadings_pc1_all_40_50[14],loadings_pc1_all_40_50[18],loadings_pc1_all_40_50[22],loadings_pc1_all_40_50[26],loadings_pc1_all_40_50[30],loadings_pc1_all_40_50[34],loadings_pc1_all_40_50[38],loadings_pc1_all_40_50[42],loadings_pc1_all_40_50[46],loadings_pc1_all_40_50[50],loadings_pc1_all_40_50[54],loadings_pc1_all_40_50[58],loadings_pc1_all_40_50[62])

loadings_pc1_all_40_50_G<-cbind(loadings_pc1_all_40_50[3],loadings_pc1_all_40_50[7],loadings_pc1_all_40_50[11],loadings_pc1_all_40_50[15],loadings_pc1_all_40_50[19],loadings_pc1_all_40_50[23],loadings_pc1_all_40_50[27],loadings_pc1_all_40_50[31],loadings_pc1_all_40_50[35],loadings_pc1_all_40_50[39],loadings_pc1_all_40_50[43],loadings_pc1_all_40_50[47],loadings_pc1_all_40_50[51],loadings_pc1_all_40_50[55],loadings_pc1_all_40_50[59],loadings_pc1_all_40_50[63])

loadings_pc1_all_40_50_U<-cbind(loadings_pc1_all_40_50[4],loadings_pc1_all_40_50[8],loadings_pc1_all_40_50[12],loadings_pc1_all_40_50[16],loadings_pc1_all_40_50[20],loadings_pc1_all_40_50[24],loadings_pc1_all_40_50[28],loadings_pc1_all_40_50[32],loadings_pc1_all_40_50[36],loadings_pc1_all_40_50[40],loadings_pc1_all_40_50[44],loadings_pc1_all_40_50[48],loadings_pc1_all_40_50[52],loadings_pc1_all_40_50[56],loadings_pc1_all_40_50[60],loadings_pc1_all_40_50[64])

# loadings pc2

loadings_pc2_all_40_50_A<-cbind(loadings_pc2_all_40_50[1],loadings_pc2_all_40_50[5],loadings_pc2_all_40_50[9],loadings_pc2_all_40_50[13],loadings_pc2_all_40_50[17],loadings_pc2_all_40_50[21],loadings_pc2_all_40_50[25],loadings_pc2_all_40_50[29],loadings_pc2_all_40_50[33],loadings_pc2_all_40_50[37],loadings_pc2_all_40_50[41],loadings_pc2_all_40_50[45],loadings_pc2_all_40_50[49],loadings_pc2_all_40_50[53],loadings_pc2_all_40_50[57],loadings_pc2_all_40_50[61])

loadings_pc2_all_40_50_C<-cbind(loadings_pc2_all_40_50[2],loadings_pc2_all_40_50[6],loadings_pc2_all_40_50[10],loadings_pc2_all_40_50[14],loadings_pc2_all_40_50[18],loadings_pc2_all_40_50[22],loadings_pc2_all_40_50[26],loadings_pc2_all_40_50[30],loadings_pc2_all_40_50[34],loadings_pc2_all_40_50[38],loadings_pc2_all_40_50[42],loadings_pc2_all_40_50[46],loadings_pc2_all_40_50[50],loadings_pc2_all_40_50[54],loadings_pc2_all_40_50[58],loadings_pc2_all_40_50[62])

loadings_pc2_all_40_50_G<-cbind(loadings_pc2_all_40_50[3],loadings_pc2_all_40_50[7],loadings_pc2_all_40_50[11],loadings_pc2_all_40_50[15],loadings_pc2_all_40_50[19],loadings_pc2_all_40_50[23],loadings_pc2_all_40_50[27],loadings_pc2_all_40_50[31],loadings_pc2_all_40_50[35],loadings_pc2_all_40_50[39],loadings_pc2_all_40_50[43],loadings_pc2_all_40_50[47],loadings_pc2_all_40_50[51],loadings_pc2_all_40_50[55],loadings_pc2_all_40_50[59],loadings_pc2_all_40_50[63])

loadings_pc2_all_40_50_U<-cbind(loadings_pc2_all_40_50[4],loadings_pc2_all_40_50[8],loadings_pc2_all_40_50[12],loadings_pc2_all_40_50[16],loadings_pc2_all_40_50[20],loadings_pc2_all_40_50[24],loadings_pc2_all_40_50[28],loadings_pc2_all_40_50[32],loadings_pc2_all_40_50[36],loadings_pc2_all_40_50[40],loadings_pc2_all_40_50[44],loadings_pc2_all_40_50[48],loadings_pc2_all_40_50[52],loadings_pc2_all_40_50[56],loadings_pc2_all_40_50[60],loadings_pc2_all_40_50[64])


text(loadings_pc1_all_40_50_A,loadings_pc2_all_40_50_A, labels_A, col="blue", xpd=T)
text(loadings_pc1_all_40_50_C,loadings_pc2_all_40_50_C, labels_C, col="red", xpd=T)
text(loadings_pc1_all_40_50_G,loadings_pc2_all_40_50_G, labels_G, col="orange", xpd=T)
text(loadings_pc1_all_40_50_U,loadings_pc2_all_40_50_U, labels_U, col="purple", xpd=T)
legend('bottomright', cex=1,ncol=1,c('A-ended','U-ended','G-ended','C-ended'),col=c('blue','purple','orange','red'),pch=c(15,15,15,15), box.col="white")
dev.off()

# e) 50-60%
pdf(file="loadings_by_kingdom_RGF_50_60.pdf", height=7, width=7)
plot(loadings_pc1_all_50_60,loadings_pc2_all_50_60, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (GC = 50-60%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)

# loadings pc1
loadings_pc1_all_50_60_A<-cbind(loadings_pc1_all_50_60[1],loadings_pc1_all_50_60[5],loadings_pc1_all_50_60[9],loadings_pc1_all_50_60[13],loadings_pc1_all_50_60[17],loadings_pc1_all_50_60[21],loadings_pc1_all_50_60[25],loadings_pc1_all_50_60[29],loadings_pc1_all_50_60[33],loadings_pc1_all_50_60[37],loadings_pc1_all_50_60[41],loadings_pc1_all_50_60[45],loadings_pc1_all_50_60[49],loadings_pc1_all_50_60[53],loadings_pc1_all_50_60[57],loadings_pc1_all_50_60[61])

loadings_pc1_all_50_60_C<-cbind(loadings_pc1_all_50_60[2],loadings_pc1_all_50_60[6],loadings_pc1_all_50_60[10],loadings_pc1_all_50_60[14],loadings_pc1_all_50_60[18],loadings_pc1_all_50_60[22],loadings_pc1_all_50_60[26],loadings_pc1_all_50_60[30],loadings_pc1_all_50_60[34],loadings_pc1_all_50_60[38],loadings_pc1_all_50_60[42],loadings_pc1_all_50_60[46],loadings_pc1_all_50_60[50],loadings_pc1_all_50_60[54],loadings_pc1_all_50_60[58],loadings_pc1_all_50_60[62])

loadings_pc1_all_50_60_G<-cbind(loadings_pc1_all_50_60[3],loadings_pc1_all_50_60[7],loadings_pc1_all_50_60[11],loadings_pc1_all_50_60[15],loadings_pc1_all_50_60[19],loadings_pc1_all_50_60[23],loadings_pc1_all_50_60[27],loadings_pc1_all_50_60[31],loadings_pc1_all_50_60[35],loadings_pc1_all_50_60[39],loadings_pc1_all_50_60[43],loadings_pc1_all_50_60[47],loadings_pc1_all_50_60[51],loadings_pc1_all_50_60[55],loadings_pc1_all_50_60[59],loadings_pc1_all_50_60[63])

loadings_pc1_all_50_60_U<-cbind(loadings_pc1_all_50_60[4],loadings_pc1_all_50_60[8],loadings_pc1_all_50_60[12],loadings_pc1_all_50_60[16],loadings_pc1_all_50_60[20],loadings_pc1_all_50_60[24],loadings_pc1_all_50_60[28],loadings_pc1_all_50_60[32],loadings_pc1_all_50_60[36],loadings_pc1_all_50_60[40],loadings_pc1_all_50_60[44],loadings_pc1_all_50_60[48],loadings_pc1_all_50_60[52],loadings_pc1_all_50_60[56],loadings_pc1_all_50_60[60],loadings_pc1_all_50_60[64])

# loadings pc2

loadings_pc2_all_50_60_A<-cbind(loadings_pc2_all_50_60[1],loadings_pc2_all_50_60[5],loadings_pc2_all_50_60[9],loadings_pc2_all_50_60[13],loadings_pc2_all_50_60[17],loadings_pc2_all_50_60[21],loadings_pc2_all_50_60[25],loadings_pc2_all_50_60[29],loadings_pc2_all_50_60[33],loadings_pc2_all_50_60[37],loadings_pc2_all_50_60[41],loadings_pc2_all_50_60[45],loadings_pc2_all_50_60[49],loadings_pc2_all_50_60[53],loadings_pc2_all_50_60[57],loadings_pc2_all_50_60[61])

loadings_pc2_all_50_60_C<-cbind(loadings_pc2_all_50_60[2],loadings_pc2_all_50_60[6],loadings_pc2_all_50_60[10],loadings_pc2_all_50_60[14],loadings_pc2_all_50_60[18],loadings_pc2_all_50_60[22],loadings_pc2_all_50_60[26],loadings_pc2_all_50_60[30],loadings_pc2_all_50_60[34],loadings_pc2_all_50_60[38],loadings_pc2_all_50_60[42],loadings_pc2_all_50_60[46],loadings_pc2_all_50_60[50],loadings_pc2_all_50_60[54],loadings_pc2_all_50_60[58],loadings_pc2_all_50_60[62])

loadings_pc2_all_50_60_G<-cbind(loadings_pc2_all_50_60[3],loadings_pc2_all_50_60[7],loadings_pc2_all_50_60[11],loadings_pc2_all_50_60[15],loadings_pc2_all_50_60[19],loadings_pc2_all_50_60[23],loadings_pc2_all_50_60[27],loadings_pc2_all_50_60[31],loadings_pc2_all_50_60[35],loadings_pc2_all_50_60[39],loadings_pc2_all_50_60[43],loadings_pc2_all_50_60[47],loadings_pc2_all_50_60[51],loadings_pc2_all_50_60[55],loadings_pc2_all_50_60[59],loadings_pc2_all_50_60[63])

loadings_pc2_all_50_60_U<-cbind(loadings_pc2_all_50_60[4],loadings_pc2_all_50_60[8],loadings_pc2_all_50_60[12],loadings_pc2_all_50_60[16],loadings_pc2_all_50_60[20],loadings_pc2_all_50_60[24],loadings_pc2_all_50_60[28],loadings_pc2_all_50_60[32],loadings_pc2_all_50_60[36],loadings_pc2_all_50_60[40],loadings_pc2_all_50_60[44],loadings_pc2_all_50_60[48],loadings_pc2_all_50_60[52],loadings_pc2_all_50_60[56],loadings_pc2_all_50_60[60],loadings_pc2_all_50_60[64])


text(loadings_pc1_all_50_60_A,loadings_pc2_all_50_60_A, labels_A, col="blue", xpd=T)
text(loadings_pc1_all_50_60_C,loadings_pc2_all_50_60_C, labels_C, col="red", xpd=T)
text(loadings_pc1_all_50_60_G,loadings_pc2_all_50_60_G, labels_G, col="orange", xpd=T)
text(loadings_pc1_all_50_60_U,loadings_pc2_all_50_60_U, labels_U, col="purple", xpd=T)
legend('bottomright', cex=1,ncol=1,c('A-ended','U-ended','G-ended','C-ended'),col=c('blue','purple','orange','red'),pch=c(15,15,15,15), box.col="white")
dev.off()


# f) 60-70%
pdf(file="loadings_by_kingdom_RGF_60_70.pdf", height=7, width=7)
plot(loadings_pc1_all_60_70,loadings_pc2_all_60_70, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (GC = 60-70%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)

# loadings pc1
loadings_pc1_all_60_70_A<-cbind(loadings_pc1_all_60_70[1],loadings_pc1_all_60_70[5],loadings_pc1_all_60_70[9],loadings_pc1_all_60_70[13],loadings_pc1_all_60_70[17],loadings_pc1_all_60_70[21],loadings_pc1_all_60_70[25],loadings_pc1_all_60_70[29],loadings_pc1_all_60_70[33],loadings_pc1_all_60_70[37],loadings_pc1_all_60_70[41],loadings_pc1_all_60_70[45],loadings_pc1_all_60_70[49],loadings_pc1_all_60_70[53],loadings_pc1_all_60_70[57],loadings_pc1_all_60_70[61])

loadings_pc1_all_60_70_C<-cbind(loadings_pc1_all_60_70[2],loadings_pc1_all_60_70[6],loadings_pc1_all_60_70[10],loadings_pc1_all_60_70[14],loadings_pc1_all_60_70[18],loadings_pc1_all_60_70[22],loadings_pc1_all_60_70[26],loadings_pc1_all_60_70[30],loadings_pc1_all_60_70[34],loadings_pc1_all_60_70[38],loadings_pc1_all_60_70[42],loadings_pc1_all_60_70[46],loadings_pc1_all_60_70[50],loadings_pc1_all_60_70[54],loadings_pc1_all_60_70[58],loadings_pc1_all_60_70[62])

loadings_pc1_all_60_70_G<-cbind(loadings_pc1_all_60_70[3],loadings_pc1_all_60_70[7],loadings_pc1_all_60_70[11],loadings_pc1_all_60_70[15],loadings_pc1_all_60_70[19],loadings_pc1_all_60_70[23],loadings_pc1_all_60_70[27],loadings_pc1_all_60_70[31],loadings_pc1_all_60_70[35],loadings_pc1_all_60_70[39],loadings_pc1_all_60_70[43],loadings_pc1_all_60_70[47],loadings_pc1_all_60_70[51],loadings_pc1_all_60_70[55],loadings_pc1_all_60_70[59],loadings_pc1_all_60_70[63])

loadings_pc1_all_60_70_U<-cbind(loadings_pc1_all_60_70[4],loadings_pc1_all_60_70[8],loadings_pc1_all_60_70[12],loadings_pc1_all_60_70[16],loadings_pc1_all_60_70[20],loadings_pc1_all_60_70[24],loadings_pc1_all_60_70[28],loadings_pc1_all_60_70[32],loadings_pc1_all_60_70[36],loadings_pc1_all_60_70[40],loadings_pc1_all_60_70[44],loadings_pc1_all_60_70[48],loadings_pc1_all_60_70[52],loadings_pc1_all_60_70[56],loadings_pc1_all_60_70[60],loadings_pc1_all_60_70[64])

# loadings pc2

loadings_pc2_all_60_70_A<-cbind(loadings_pc2_all_60_70[1],loadings_pc2_all_60_70[5],loadings_pc2_all_60_70[9],loadings_pc2_all_60_70[13],loadings_pc2_all_60_70[17],loadings_pc2_all_60_70[21],loadings_pc2_all_60_70[25],loadings_pc2_all_60_70[29],loadings_pc2_all_60_70[33],loadings_pc2_all_60_70[37],loadings_pc2_all_60_70[41],loadings_pc2_all_60_70[45],loadings_pc2_all_60_70[49],loadings_pc2_all_60_70[53],loadings_pc2_all_60_70[57],loadings_pc2_all_60_70[61])

loadings_pc2_all_60_70_C<-cbind(loadings_pc2_all_60_70[2],loadings_pc2_all_60_70[6],loadings_pc2_all_60_70[10],loadings_pc2_all_60_70[14],loadings_pc2_all_60_70[18],loadings_pc2_all_60_70[22],loadings_pc2_all_60_70[26],loadings_pc2_all_60_70[30],loadings_pc2_all_60_70[34],loadings_pc2_all_60_70[38],loadings_pc2_all_60_70[42],loadings_pc2_all_60_70[46],loadings_pc2_all_60_70[50],loadings_pc2_all_60_70[54],loadings_pc2_all_60_70[58],loadings_pc2_all_60_70[62])

loadings_pc2_all_60_70_G<-cbind(loadings_pc2_all_60_70[3],loadings_pc2_all_60_70[7],loadings_pc2_all_60_70[11],loadings_pc2_all_60_70[15],loadings_pc2_all_60_70[19],loadings_pc2_all_60_70[23],loadings_pc2_all_60_70[27],loadings_pc2_all_60_70[31],loadings_pc2_all_60_70[35],loadings_pc2_all_60_70[39],loadings_pc2_all_60_70[43],loadings_pc2_all_60_70[47],loadings_pc2_all_60_70[51],loadings_pc2_all_60_70[55],loadings_pc2_all_60_70[59],loadings_pc2_all_60_70[63])

loadings_pc2_all_60_70_U<-cbind(loadings_pc2_all_60_70[4],loadings_pc2_all_60_70[8],loadings_pc2_all_60_70[12],loadings_pc2_all_60_70[16],loadings_pc2_all_60_70[20],loadings_pc2_all_60_70[24],loadings_pc2_all_60_70[28],loadings_pc2_all_60_70[32],loadings_pc2_all_60_70[36],loadings_pc2_all_60_70[40],loadings_pc2_all_60_70[44],loadings_pc2_all_60_70[48],loadings_pc2_all_60_70[52],loadings_pc2_all_60_70[56],loadings_pc2_all_60_70[60],loadings_pc2_all_60_70[64])


text(loadings_pc1_all_60_70_A,loadings_pc2_all_60_70_A, labels_A, col="blue", xpd=T)
text(loadings_pc1_all_60_70_C,loadings_pc2_all_60_70_C, labels_C, col="red", xpd=T)
text(loadings_pc1_all_60_70_G,loadings_pc2_all_60_70_G, labels_G, col="orange", xpd=T)
text(loadings_pc1_all_60_70_U,loadings_pc2_all_60_70_U, labels_U, col="purple", xpd=T)
legend('bottomright', cex=1,ncol=1,c('A-ended','U-ended','G-ended','C-ended'),col=c('blue','purple','orange','red'),pch=c(15,15,15,15), box.col="white")
dev.off()



##########################

# ANALYSIS OF ONLY ARG ###

##########################

arg_dat_allGC<-cbind(data_allGC$CGU, data_allGC$CGC, data_allGC$CGA, data_allGC$CGG, data_allGC$AGA, data_allGC$AGG)
arg_dat_below30<-cbind(data_below30$CGU, data_below30$CGC, data_below30$CGA, data_below30$CGG, data_below30$AGA, data_below30$AGG)
arg_dat_30_40<-cbind(data_30_40$CGU, data_30_40$CGC, data_30_40$CGA, data_30_40$CGG, data_30_40$AGA, data_30_40$AGG)
arg_dat_40_50<-cbind(data_40_50$CGU, data_40_50$CGC, data_40_50$CGA, data_40_50$CGG, data_40_50$AGA, data_40_50$AGG)
arg_dat_50_60<-cbind(data_50_60$CGU, data_50_60$CGC, data_50_60$CGA, data_50_60$CGG, data_50_60$AGA, data_50_60$AGG)
arg_dat_60_70<-cbind(data_60_70$CGU, data_60_70$CGC, data_60_70$CGA, data_60_70$CGG, data_60_70$AGA, data_60_70$AGG)

colnames(arg_dat_allGC)=c("CGU","CGC","CGA","CGG","AGA","AGG")
colnames(arg_dat_below30)=c("CGU","CGC","CGA","CGG","AGA","AGG")
colnames(arg_dat_30_40)=c("CGU","CGC","CGA","CGG","AGA","AGG")
colnames(arg_dat_40_50)=c("CGU","CGC","CGA","CGG","AGA","AGG")
colnames(arg_dat_50_60)=c("CGU","CGC","CGA","CGG","AGA","AGG")
colnames(arg_dat_60_70)=c("CGU","CGC","CGA","CGG","AGA","AGG")

arg_model_allGC<-prcomp(arg_dat_allGC)
arg_model_below30<-prcomp(arg_dat_below30)
arg_model_30_40<-prcomp(arg_dat_30_40)
arg_model_40_50<-prcomp(arg_dat_40_50)
arg_model_50_60<-prcomp(arg_dat_50_60)
arg_model_60_70<-prcomp(arg_dat_60_70)

arg_scores_allGC<-arg_model_allGC$x
arg_scores_below30<-arg_model_below30$x
arg_scores_30_40<-arg_model_30_40$x
arg_scores_40_50<-arg_model_40_50$x
arg_scores_50_60<-arg_model_50_60$x
arg_scores_60_70<-arg_model_60_70$x

arg_scores_data_allGC<-as.matrix(arg_scores_allGC[,1:2])
arg_scores_data_below30<-as.matrix(arg_scores_below30[,1:2])
arg_scores_data_30_40<-as.matrix(arg_scores_30_40[,1:2])
arg_scores_data_40_50<-as.matrix(arg_scores_40_50[,1:2])
arg_scores_data_50_60<-as.matrix(arg_scores_50_60[,1:2])
arg_scores_data_60_70<-as.matrix(arg_scores_60_70[,1:2])

arg_pc1_all_allGC<-as.vector(arg_scores_data_allGC[,1])
arg_pc2_all_allGC<-as.vector(arg_scores_data_allGC[,2])
arg_pc1_all_below30<-as.vector(arg_scores_data_below30[,1])
arg_pc2_all_below30<-as.vector(arg_scores_data_below30[,2])
arg_pc1_all_30_40<-as.vector(arg_scores_data_30_40[,1])
arg_pc2_all_30_40<-as.vector(arg_scores_data_30_40[,2])
arg_pc1_all_40_50<-as.vector(arg_scores_data_40_50[,1])
arg_pc2_all_40_50<-as.vector(arg_scores_data_40_50[,2])
arg_pc1_all_50_60<-as.vector(arg_scores_data_50_60[,1])
arg_pc2_all_50_60<-as.vector(arg_scores_data_50_60[,2])
arg_pc1_all_60_70<-as.vector(arg_scores_data_60_70[,1])
arg_pc2_all_60_70<-as.vector(arg_scores_data_60_70[,2])

#GET SCORES FOR EUK AND BACT

arg_arch_allGC<-arg_scores_data_allGC[1:120,] 		# 120
arg_bact_allGC<-arg_scores_data_allGC[121:1461,]  	# 1341
arg_euk_allGC<-arg_scores_data_allGC[1462:1625,]	# 164 --> total = 1625
arg_pc1_arch_allGC<-as.vector(arg_arch_allGC[,1])
arg_pc2_arch_allGC<-as.vector(arg_arch_allGC[,2])
arg_pc1_euk_allGC<-as.vector(arg_euk_allGC[,1])
arg_pc2_euk_allGC<-as.vector(arg_euk_allGC[,2])
arg_pc1_bact_allGC<-as.vector(arg_bact_allGC[,1])
arg_pc2_bact_allGC<-as.vector(arg_bact_allGC[,2])

#arg_arch_below30<-arg_scores_data_below30[1:1,] 	# 1
arg_bact_below30<-arg_scores_data_below30[2:37,]  	# 36
arg_euk_below30<-arg_scores_data_below30[38:47,]	# 10 --> total = 47
#arg_pc1_arch_below30<-as.vector(arg_arch_below30[,1])
#arg_pc2_arch_below30<-as.vector(arg_arch_below30[,2])
arg_pc1_euk_below30<-as.vector(arg_euk_below30[,1])
arg_pc2_euk_below30<-as.vector(arg_euk_below30[,2])
arg_pc1_bact_below30<-as.vector(arg_bact_below30[,1])
arg_pc2_bact_below30<-as.vector(arg_bact_below30[,2])

arg_arch_30_40<-arg_scores_data_30_40[1:26,]	# 26
arg_bact_30_40<-arg_scores_data_30_40[27:331,]  # 305
arg_euk_30_40<-arg_scores_data_30_40[332:356,]	# 25 --> total = 356
arg_pc1_arch_30_40<-as.vector(arg_arch_30_40[,1])
arg_pc2_arch_30_40<-as.vector(arg_arch_30_40[,2])
arg_pc1_euk_30_40<-as.vector(arg_euk_30_40[,1])
arg_pc2_euk_30_40<-as.vector(arg_euk_30_40[,2])
arg_pc1_bact_30_40<-as.vector(arg_bact_30_40[,1])
arg_pc2_bact_30_40<-as.vector(arg_bact_30_40[,2])

arg_arch_40_50<-arg_scores_data_40_50[1:35,]	# 35
arg_bact_40_50<-arg_scores_data_40_50[36:375,]  # 340
arg_euk_40_50<-arg_scores_data_40_50[376:439,]	# 64 --> total = 439
arg_pc1_arch_40_50<-as.vector(arg_arch_40_50[,1])
arg_pc2_arch_40_50<-as.vector(arg_arch_40_50[,2])
arg_pc1_euk_40_50<-as.vector(arg_euk_40_50[,1])
arg_pc2_euk_40_50<-as.vector(arg_euk_40_50[,2])
arg_pc1_bact_40_50<-as.vector(arg_bact_40_50[,1])
arg_pc2_bact_40_50<-as.vector(arg_bact_40_50[,2])

arg_arch_50_60<-arg_scores_data_50_60[1:34,]	# 34
arg_bact_50_60<-arg_scores_data_50_60[35:327,]  # 293
arg_euk_50_60<-arg_scores_data_50_60[328:377,]	# 50 --> total = 377
arg_pc1_arch_50_60<-as.vector(arg_arch_50_60[,1])
arg_pc2_arch_50_60<-as.vector(arg_arch_50_60[,2])
arg_pc1_euk_50_60<-as.vector(arg_euk_50_60[,1])
arg_pc2_euk_50_60<-as.vector(arg_euk_50_60[,2])
arg_pc1_bact_50_60<-as.vector(arg_bact_50_60[,1])
arg_pc2_bact_50_60<-as.vector(arg_bact_50_60[,2])

arg_arch_60_70<-arg_scores_data_60_70[1:24,]	# 24
arg_bact_60_70<-arg_scores_data_60_70[25:326,]  # 302
arg_euk_60_70<-arg_scores_data_60_70[327:340,]	# 14 --> total = 340
arg_pc1_arch_60_70<-as.vector(arg_arch_60_70[,1])
arg_pc2_arch_60_70<-as.vector(arg_arch_60_70[,2])
arg_pc1_euk_60_70<-as.vector(arg_euk_60_70[,1])
arg_pc2_euk_60_70<-as.vector(arg_euk_60_70[,2])
arg_pc1_bact_60_70<-as.vector(arg_bact_60_70[,1])
arg_pc2_bact_60_70<-as.vector(arg_bact_60_70[,2])

# PLOT PCA SCORES

pdf(file="arg_scores_by_kingdom_RGF_allGC.pdf", height=7, width=7)
plot(arg_pc1_all_allGC, arg_pc2_all_allGC, type="n", xlab="PC1", ylab="PC2")
title(main="ARG codon usage across kingdoms (all GC's)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(arg_pc1_arch_allGC, arg_pc2_arch_allGC, pch=15, col="red")
points(arg_pc1_bact_allGC, arg_pc2_bact_allGC, pch=16, col="purple")
points(arg_pc1_euk_allGC, arg_pc2_euk_allGC, pch=17,col="forestgreen")
#arrows(0, 0, X[,1], X[,2], len=0.1, col="black")
#text(1.4*X, rownames(X), col="black", xpd=T)
legend('topright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="arg_scores_by_kingdom_RGF_below30.pdf", height=7, width=7)
plot(arg_pc1_all_below30, arg_pc2_all_below30, type="n", xlab="PC1", ylab="PC2")
title(main="ARG codon usage across kingdoms (GC < 30%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
#points(arg_pc1_arch_below30, arg_pc2_arch_below30, pch=15, col="red")
points(arg_pc1_bact_below30, arg_pc2_bact_below30, pch=16, col="purple")
points(arg_pc1_euk_below30, arg_pc2_euk_below30, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Bacteria','Eukarya'),col=c('purple','forestgreen'),pch=c(16,17), box.col="white")
dev.off()

pdf(file="arg_scores_by_kingdom_RGF_30_40.pdf", height=7, width=7)
plot(arg_pc1_all_30_40, arg_pc2_all_30_40, type="n", xlab="PC1", ylab="PC2")
title(main="ARG codon usage across kingdoms (GC = 30-40%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(arg_pc1_arch_30_40, arg_pc2_arch_30_40, pch=15, col="red")
points(arg_pc1_bact_30_40, arg_pc2_bact_30_40, pch=16, col="purple")
points(arg_pc1_euk_30_40, arg_pc2_euk_30_40, pch=17,col="forestgreen")
legend('bottomright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="arg_scores_by_kingdom_RGF_40_50.pdf", height=7, width=7)
plot(arg_pc1_all_40_50, arg_pc2_all_40_50, type="n", xlab="PC1", ylab="PC2")
title(main="ARG codon usage across kingdoms (GC = 40-50%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(arg_pc1_arch_40_50, arg_pc2_arch_40_50, pch=15, col="red")
points(arg_pc1_bact_40_50, arg_pc2_bact_40_50, pch=16, col="purple")
points(arg_pc1_euk_40_50, arg_pc2_euk_40_50, pch=17,col="forestgreen")
legend('bottomright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="arg_scores_by_kingdom_RGF_50_60.pdf", height=7, width=7)
plot(arg_pc1_all_50_60, arg_pc2_all_50_60, type="n", xlab="PC1", ylab="PC2")
title(main="ARG codon usage across kingdoms (GC = 50-60%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(arg_pc1_arch_50_60, arg_pc2_arch_50_60, pch=15, col="red")
points(arg_pc1_bact_50_60, arg_pc2_bact_50_60, pch=16, col="purple")
points(arg_pc1_euk_50_60, arg_pc2_euk_50_60, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="arg_scores_by_kingdom_RGF_60_70.pdf", height=7, width=7)
plot(arg_pc1_all_60_70, arg_pc2_all_60_70, type="n", xlab="PC1", ylab="PC2")
title(main="ARG codon usage across kingdoms (GC = 60-70%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(arg_pc1_arch_60_70, arg_pc2_arch_60_70, pch=15, col="red")
points(arg_pc1_bact_60_70, arg_pc2_bact_60_70, pch=16, col="purple")
points(arg_pc1_euk_60_70, arg_pc2_euk_60_70, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

##########################

# ANALYSIS OF ONLY LEU ###

##########################

leu_dat_allGC<-cbind(data_allGC$CUU, data_allGC$CUC, data_allGC$CUA, data_allGC$CUG, data_allGC$UUA, data_allGC$UUG)
leu_dat_below30<-cbind(data_below30$CUU, data_below30$CUC, data_below30$CUA, data_below30$CUG, data_below30$UUA, data_below30$UUG)
leu_dat_30_40<-cbind(data_30_40$CUU, data_30_40$CUC, data_30_40$CUA, data_30_40$CUG, data_30_40$UUA, data_30_40$UUG)
leu_dat_40_50<-cbind(data_40_50$CUU, data_40_50$CUC, data_40_50$CUA, data_40_50$CUG, data_40_50$UUA, data_40_50$UUG)
leu_dat_50_60<-cbind(data_50_60$CUU, data_50_60$CUC, data_50_60$CUA, data_50_60$CUG, data_50_60$UUA, data_50_60$UUG)
leu_dat_60_70<-cbind(data_60_70$CUU, data_60_70$CUC, data_60_70$CUA, data_60_70$CUG, data_60_70$UUA, data_60_70$UUG)

colnames(leu_dat_allGC)=c("CUU","CUC","CUA","CUG","UUA","UUG")
colnames(leu_dat_below30)=c("CUU","CUC","CUA","CUG","UUA","UUG")
colnames(leu_dat_30_40)=c("CUU","CUC","CUA","CUG","UUA","UUG")
colnames(leu_dat_40_50)=c("CUU","CUC","CUA","CUG","UUA","UUG")
colnames(leu_dat_50_60)=c("CUU","CUC","CUA","CUG","UUA","UUG")
colnames(leu_dat_60_70)=c("CUU","CUC","CUA","CUG","UUA","UUG")

leu_model_allGC<-prcomp(leu_dat_allGC)
leu_model_below30<-prcomp(leu_dat_below30)
leu_model_30_40<-prcomp(leu_dat_30_40)
leu_model_40_50<-prcomp(leu_dat_40_50)
leu_model_50_60<-prcomp(leu_dat_50_60)
leu_model_60_70<-prcomp(leu_dat_60_70)

leu_scores_allGC<-leu_model_allGC$x
leu_scores_below30<-leu_model_below30$x
leu_scores_30_40<-leu_model_30_40$x
leu_scores_40_50<-leu_model_40_50$x
leu_scores_50_60<-leu_model_50_60$x
leu_scores_60_70<-leu_model_60_70$x

leu_scores_data_allGC<-as.matrix(leu_scores_allGC[,1:2])
leu_scores_data_below30<-as.matrix(leu_scores_below30[,1:2])
leu_scores_data_30_40<-as.matrix(leu_scores_30_40[,1:2])
leu_scores_data_40_50<-as.matrix(leu_scores_40_50[,1:2])
leu_scores_data_50_60<-as.matrix(leu_scores_50_60[,1:2])
leu_scores_data_60_70<-as.matrix(leu_scores_60_70[,1:2])

leu_pc1_all_allGC<-as.vector(leu_scores_data_allGC[,1])
leu_pc2_all_allGC<-as.vector(leu_scores_data_allGC[,2])
leu_pc1_all_below30<-as.vector(leu_scores_data_below30[,1])
leu_pc2_all_below30<-as.vector(leu_scores_data_below30[,2])
leu_pc1_all_30_40<-as.vector(leu_scores_data_30_40[,1])
leu_pc2_all_30_40<-as.vector(leu_scores_data_30_40[,2])
leu_pc1_all_40_50<-as.vector(leu_scores_data_40_50[,1])
leu_pc2_all_40_50<-as.vector(leu_scores_data_40_50[,2])
leu_pc1_all_50_60<-as.vector(leu_scores_data_50_60[,1])
leu_pc2_all_50_60<-as.vector(leu_scores_data_50_60[,2])
leu_pc1_all_60_70<-as.vector(leu_scores_data_60_70[,1])
leu_pc2_all_60_70<-as.vector(leu_scores_data_60_70[,2])

#GET SCORES FOR EUK AND BACT

leu_arch_allGC<-leu_scores_data_allGC[1:120,] 		# 120
leu_bact_allGC<-leu_scores_data_allGC[121:1461,]  	# 1341
leu_euk_allGC<-leu_scores_data_allGC[1462:1625,]	# 164 --> total = 1625
leu_pc1_arch_allGC<-as.vector(leu_arch_allGC[,1])
leu_pc2_arch_allGC<-as.vector(leu_arch_allGC[,2])
leu_pc1_euk_allGC<-as.vector(leu_euk_allGC[,1])
leu_pc2_euk_allGC<-as.vector(leu_euk_allGC[,2])
leu_pc1_bact_allGC<-as.vector(leu_bact_allGC[,1])
leu_pc2_bact_allGC<-as.vector(leu_bact_allGC[,2])

#leu_arch_below30<-leu_scores_data_below30[1:1,] 	# 1
leu_bact_below30<-leu_scores_data_below30[2:37,]  	# 36
leu_euk_below30<-leu_scores_data_below30[38:47,]	# 10 --> total = 47
#leu_pc1_arch_below30<-as.vector(leu_arch_below30[,1])
#leu_pc2_arch_below30<-as.vector(leu_arch_below30[,2])
leu_pc1_euk_below30<-as.vector(leu_euk_below30[,1])
leu_pc2_euk_below30<-as.vector(leu_euk_below30[,2])
leu_pc1_bact_below30<-as.vector(leu_bact_below30[,1])
leu_pc2_bact_below30<-as.vector(leu_bact_below30[,2])

leu_arch_30_40<-leu_scores_data_30_40[1:26,]	# 26
leu_bact_30_40<-leu_scores_data_30_40[27:331,]  # 305
leu_euk_30_40<-leu_scores_data_30_40[332:356,]	# 25 --> total = 356
leu_pc1_arch_30_40<-as.vector(leu_arch_30_40[,1])
leu_pc2_arch_30_40<-as.vector(leu_arch_30_40[,2])
leu_pc1_euk_30_40<-as.vector(leu_euk_30_40[,1])
leu_pc2_euk_30_40<-as.vector(leu_euk_30_40[,2])
leu_pc1_bact_30_40<-as.vector(leu_bact_30_40[,1])
leu_pc2_bact_30_40<-as.vector(leu_bact_30_40[,2])

leu_arch_40_50<-leu_scores_data_40_50[1:35,]	# 35
leu_bact_40_50<-leu_scores_data_40_50[36:375,]  # 340
leu_euk_40_50<-leu_scores_data_40_50[376:439,]	# 64 --> total = 439
leu_pc1_arch_40_50<-as.vector(leu_arch_40_50[,1])
leu_pc2_arch_40_50<-as.vector(leu_arch_40_50[,2])
leu_pc1_euk_40_50<-as.vector(leu_euk_40_50[,1])
leu_pc2_euk_40_50<-as.vector(leu_euk_40_50[,2])
leu_pc1_bact_40_50<-as.vector(leu_bact_40_50[,1])
leu_pc2_bact_40_50<-as.vector(leu_bact_40_50[,2])

leu_arch_50_60<-leu_scores_data_50_60[1:34,]	# 34
leu_bact_50_60<-leu_scores_data_50_60[35:327,]  # 293
leu_euk_50_60<-leu_scores_data_50_60[328:377,]	# 50 --> total = 377
leu_pc1_arch_50_60<-as.vector(leu_arch_50_60[,1])
leu_pc2_arch_50_60<-as.vector(leu_arch_50_60[,2])
leu_pc1_euk_50_60<-as.vector(leu_euk_50_60[,1])
leu_pc2_euk_50_60<-as.vector(leu_euk_50_60[,2])
leu_pc1_bact_50_60<-as.vector(leu_bact_50_60[,1])
leu_pc2_bact_50_60<-as.vector(leu_bact_50_60[,2])

leu_arch_60_70<-leu_scores_data_60_70[1:24,]	# 24
leu_bact_60_70<-leu_scores_data_60_70[25:326,]  # 302
leu_euk_60_70<-leu_scores_data_60_70[327:340,]	# 14 --> total = 340
leu_pc1_arch_60_70<-as.vector(leu_arch_60_70[,1])
leu_pc2_arch_60_70<-as.vector(leu_arch_60_70[,2])
leu_pc1_euk_60_70<-as.vector(leu_euk_60_70[,1])
leu_pc2_euk_60_70<-as.vector(leu_euk_60_70[,2])
leu_pc1_bact_60_70<-as.vector(leu_bact_60_70[,1])
leu_pc2_bact_60_70<-as.vector(leu_bact_60_70[,2])

# PLOT PCA SCORES

pdf(file="leu_scores_by_kingdom_RGF_allGC.pdf", height=7, width=7)
plot(leu_pc1_all_allGC, leu_pc2_all_allGC, type="n", xlab="PC1", ylab="PC2")
title(main="LEU codon usage across kingdoms (all GC's)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(leu_pc1_arch_allGC, leu_pc2_arch_allGC, pch=15, col="red")
points(leu_pc1_bact_allGC, leu_pc2_bact_allGC, pch=16, col="purple")
points(leu_pc1_euk_allGC, leu_pc2_euk_allGC, pch=17,col="forestgreen")
#arrows(0, 0, X[,1], X[,2], len=0.1, col="black")
#text(1.4*X, rownames(X), col="black", xpd=T)
legend('topright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="leu_scores_by_kingdom_RGF_below30.pdf", height=7, width=7)
plot(leu_pc1_all_below30, leu_pc2_all_below30, type="n", xlab="PC1", ylab="PC2")
title(main="LEU codon usage across kingdoms (GC < 30%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
#points(leu_pc1_arch_below30, leu_pc2_arch_below30, pch=15, col="red")
points(leu_pc1_bact_below30, leu_pc2_bact_below30, pch=16, col="purple")
points(leu_pc1_euk_below30, leu_pc2_euk_below30, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Bacteria','Eukarya'),col=c('purple','forestgreen'),pch=c(16,17), box.col="white")
dev.off()

pdf(file="leu_scores_by_kingdom_RGF_30_40.pdf", height=7, width=7)
plot(leu_pc1_all_30_40, leu_pc2_all_30_40, type="n", xlab="PC1", ylab="PC2")
title(main="LEU codon usage across kingdoms (GC = 30-40%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(leu_pc1_arch_30_40, leu_pc2_arch_30_40, pch=15, col="red")
points(leu_pc1_bact_30_40, leu_pc2_bact_30_40, pch=16, col="purple")
points(leu_pc1_euk_30_40, leu_pc2_euk_30_40, pch=17,col="forestgreen")
legend('bottomright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="leu_scores_by_kingdom_RGF_40_50.pdf", height=7, width=7)
plot(leu_pc1_all_40_50, leu_pc2_all_40_50, type="n", xlab="PC1", ylab="PC2")
title(main="LEU codon usage across kingdoms (GC = 40-50%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(leu_pc1_arch_40_50, leu_pc2_arch_40_50, pch=15, col="red")
points(leu_pc1_bact_40_50, leu_pc2_bact_40_50, pch=16, col="purple")
points(leu_pc1_euk_40_50, leu_pc2_euk_40_50, pch=17,col="forestgreen")
legend('bottomright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="leu_scores_by_kingdom_RGF_50_60.pdf", height=7, width=7)
plot(leu_pc1_all_50_60, leu_pc2_all_50_60, type="n", xlab="PC1", ylab="PC2")
title(main="LEU codon usage across kingdoms (GC = 50-60%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(leu_pc1_arch_50_60, leu_pc2_arch_50_60, pch=15, col="red")
points(leu_pc1_bact_50_60, leu_pc2_bact_50_60, pch=16, col="purple")
points(leu_pc1_euk_50_60, leu_pc2_euk_50_60, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="leu_scores_by_kingdom_RGF_60_70.pdf", height=7, width=7)
plot(leu_pc1_all_60_70, leu_pc2_all_60_70, type="n", xlab="PC1", ylab="PC2")
title(main="LEU codon usage across kingdoms (GC = 60-70%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(leu_pc1_arch_60_70, leu_pc2_arch_60_70, pch=15, col="red")
points(leu_pc1_bact_60_70, leu_pc2_bact_60_70, pch=16, col="purple")
points(leu_pc1_euk_60_70, leu_pc2_euk_60_70, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()


##########################

# ANALYSIS OF ONLY ARG-LEU ###

##########################

argleu_dat_allGC<-cbind(data_allGC$CGU, data_allGC$CGC, data_allGC$CGA, data_allGC$CGG, data_allGC$AGA, data_allGC$AGG,data_allGC$CUU, data_allGC$CUC, data_allGC$CUA, data_allGC$CUG, data_allGC$UUA, data_allGC$UUG)
argleu_dat_below30<-cbind(data_below30$CGU, data_below30$CGC, data_below30$CGA, data_below30$CGG, data_below30$AGA, data_below30$AGG, data_below30$CUU, data_below30$CUC, data_below30$CUA, data_below30$CUG, data_below30$UUA, data_below30$UUG)
argleu_dat_30_40<-cbind(data_30_40$CGU, data_30_40$CGC, data_30_40$CGA, data_30_40$CGG, data_30_40$AGA, data_30_40$AGG, data_30_40$CUU, data_30_40$CUC, data_30_40$CUA, data_30_40$CUG, data_30_40$UUA, data_30_40$UUG)
argleu_dat_40_50<-cbind(data_40_50$CGU, data_40_50$CGC, data_40_50$CGA, data_40_50$CGG, data_40_50$AGA, data_40_50$AGG, data_40_50$CUU, data_40_50$CUC, data_40_50$CUA, data_40_50$CUG, data_40_50$UUA, data_40_50$UUG)
argleu_dat_50_60<-cbind(data_50_60$CGU, data_50_60$CGC, data_50_60$CGA, data_50_60$CGG, data_50_60$AGA, data_50_60$AGG, data_50_60$CUU, data_50_60$CUC, data_50_60$CUA, data_50_60$CUG, data_50_60$UUA, data_50_60$UUG)
argleu_dat_60_70<-cbind(data_60_70$CGU, data_60_70$CGC, data_60_70$CGA, data_60_70$CGG, data_60_70$AGA, data_60_70$AGG, data_60_70$CUU, data_60_70$CUC, data_60_70$CUA, data_60_70$CUG, data_60_70$UUA, data_60_70$UUG)

colnames(argleu_dat_allGC)=c("CGU","CGC","CGA","CGG","AGA","AGG","CUU","CUC","CUA","CUG","UUA","UUG")
colnames(argleu_dat_below30)=c("CGU","CGC","CGA","CGG","AGA","AGG","CUU","CUC","CUA","CUG","UUA","UUG")
colnames(argleu_dat_30_40)=c("CGU","CGC","CGA","CGG","AGA","AGG","CUU","CUC","CUA","CUG","UUA","UUG")
colnames(argleu_dat_40_50)=c("CGU","CGC","CGA","CGG","AGA","AGG","CUU","CUC","CUA","CUG","UUA","UUG")
colnames(argleu_dat_50_60)=c("CGU","CGC","CGA","CGG","AGA","AGG","CUU","CUC","CUA","CUG","UUA","UUG")
colnames(argleu_dat_60_70)=c("CGU","CGC","CGA","CGG","AGA","AGG","CUU","CUC","CUA","CUG","UUA","UUG")

argleu_model_allGC<-prcomp(argleu_dat_allGC)
argleu_model_below30<-prcomp(argleu_dat_below30)
argleu_model_30_40<-prcomp(argleu_dat_30_40)
argleu_model_40_50<-prcomp(argleu_dat_40_50)
argleu_model_50_60<-prcomp(argleu_dat_50_60)
argleu_model_60_70<-prcomp(argleu_dat_60_70)

argleu_scores_allGC<-argleu_model_allGC$x
argleu_scores_below30<-argleu_model_below30$x
argleu_scores_30_40<-argleu_model_30_40$x
argleu_scores_40_50<-argleu_model_40_50$x
argleu_scores_50_60<-argleu_model_50_60$x
argleu_scores_60_70<-argleu_model_60_70$x

argleu_scores_data_allGC<-as.matrix(argleu_scores_allGC[,1:2])
argleu_scores_data_below30<-as.matrix(argleu_scores_below30[,1:2])
argleu_scores_data_30_40<-as.matrix(argleu_scores_30_40[,1:2])
argleu_scores_data_40_50<-as.matrix(argleu_scores_40_50[,1:2])
argleu_scores_data_50_60<-as.matrix(argleu_scores_50_60[,1:2])
argleu_scores_data_60_70<-as.matrix(argleu_scores_60_70[,1:2])

argleu_pc1_all_allGC<-as.vector(argleu_scores_data_allGC[,1])
argleu_pc2_all_allGC<-as.vector(argleu_scores_data_allGC[,2])
argleu_pc1_all_below30<-as.vector(argleu_scores_data_below30[,1])
argleu_pc2_all_below30<-as.vector(argleu_scores_data_below30[,2])
argleu_pc1_all_30_40<-as.vector(argleu_scores_data_30_40[,1])
argleu_pc2_all_30_40<-as.vector(argleu_scores_data_30_40[,2])
argleu_pc1_all_40_50<-as.vector(argleu_scores_data_40_50[,1])
argleu_pc2_all_40_50<-as.vector(argleu_scores_data_40_50[,2])
argleu_pc1_all_50_60<-as.vector(argleu_scores_data_50_60[,1])
argleu_pc2_all_50_60<-as.vector(argleu_scores_data_50_60[,2])
argleu_pc1_all_60_70<-as.vector(argleu_scores_data_60_70[,1])
argleu_pc2_all_60_70<-as.vector(argleu_scores_data_60_70[,2])

#GET SCORES FOR EUK AND BACT

argleu_arch_allGC<-argleu_scores_data_allGC[1:120,] 		# 120
argleu_bact_allGC<-argleu_scores_data_allGC[121:1461,]  	# 1341
argleu_euk_allGC<-argleu_scores_data_allGC[1462:1625,]	# 164 --> total = 1625
argleu_pc1_arch_allGC<-as.vector(argleu_arch_allGC[,1])
argleu_pc2_arch_allGC<-as.vector(argleu_arch_allGC[,2])
argleu_pc1_euk_allGC<-as.vector(argleu_euk_allGC[,1])
argleu_pc2_euk_allGC<-as.vector(argleu_euk_allGC[,2])
argleu_pc1_bact_allGC<-as.vector(argleu_bact_allGC[,1])
argleu_pc2_bact_allGC<-as.vector(argleu_bact_allGC[,2])

#argleu_arch_below30<-argleu_scores_data_below30[1:1,] 	# 1
argleu_bact_below30<-argleu_scores_data_below30[2:37,]  	# 36
argleu_euk_below30<-argleu_scores_data_below30[38:47,]	# 10 --> total = 47
#argleu_pc1_arch_below30<-as.vector(argleu_arch_below30[,1])
#argleu_pc2_arch_below30<-as.vector(argleu_arch_below30[,2])
argleu_pc1_euk_below30<-as.vector(argleu_euk_below30[,1])
argleu_pc2_euk_below30<-as.vector(argleu_euk_below30[,2])
argleu_pc1_bact_below30<-as.vector(argleu_bact_below30[,1])
argleu_pc2_bact_below30<-as.vector(argleu_bact_below30[,2])

argleu_arch_30_40<-argleu_scores_data_30_40[1:26,]	# 26
argleu_bact_30_40<-argleu_scores_data_30_40[27:331,]  # 305
argleu_euk_30_40<-argleu_scores_data_30_40[332:356,]	# 25 --> total = 356
argleu_pc1_arch_30_40<-as.vector(argleu_arch_30_40[,1])
argleu_pc2_arch_30_40<-as.vector(argleu_arch_30_40[,2])
argleu_pc1_euk_30_40<-as.vector(argleu_euk_30_40[,1])
argleu_pc2_euk_30_40<-as.vector(argleu_euk_30_40[,2])
argleu_pc1_bact_30_40<-as.vector(argleu_bact_30_40[,1])
argleu_pc2_bact_30_40<-as.vector(argleu_bact_30_40[,2])

argleu_arch_40_50<-argleu_scores_data_40_50[1:35,]	# 35
argleu_bact_40_50<-argleu_scores_data_40_50[36:375,]  # 340
argleu_euk_40_50<-argleu_scores_data_40_50[376:439,]	# 64 --> total = 439
argleu_pc1_arch_40_50<-as.vector(argleu_arch_40_50[,1])
argleu_pc2_arch_40_50<-as.vector(argleu_arch_40_50[,2])
argleu_pc1_euk_40_50<-as.vector(argleu_euk_40_50[,1])
argleu_pc2_euk_40_50<-as.vector(argleu_euk_40_50[,2])
argleu_pc1_bact_40_50<-as.vector(argleu_bact_40_50[,1])
argleu_pc2_bact_40_50<-as.vector(argleu_bact_40_50[,2])

argleu_arch_50_60<-argleu_scores_data_50_60[1:34,]	# 34
argleu_bact_50_60<-argleu_scores_data_50_60[35:327,]  # 293
argleu_euk_50_60<-argleu_scores_data_50_60[328:377,]	# 50 --> total = 377
argleu_pc1_arch_50_60<-as.vector(argleu_arch_50_60[,1])
argleu_pc2_arch_50_60<-as.vector(argleu_arch_50_60[,2])
argleu_pc1_euk_50_60<-as.vector(argleu_euk_50_60[,1])
argleu_pc2_euk_50_60<-as.vector(argleu_euk_50_60[,2])
argleu_pc1_bact_50_60<-as.vector(argleu_bact_50_60[,1])
argleu_pc2_bact_50_60<-as.vector(argleu_bact_50_60[,2])

argleu_arch_60_70<-argleu_scores_data_60_70[1:24,]	# 24
argleu_bact_60_70<-argleu_scores_data_60_70[25:326,]  # 302
argleu_euk_60_70<-argleu_scores_data_60_70[327:340,]	# 14 --> total = 340
argleu_pc1_arch_60_70<-as.vector(argleu_arch_60_70[,1])
argleu_pc2_arch_60_70<-as.vector(argleu_arch_60_70[,2])
argleu_pc1_euk_60_70<-as.vector(argleu_euk_60_70[,1])
argleu_pc2_euk_60_70<-as.vector(argleu_euk_60_70[,2])
argleu_pc1_bact_60_70<-as.vector(argleu_bact_60_70[,1])
argleu_pc2_bact_60_70<-as.vector(argleu_bact_60_70[,2])

# PLOT PCA SCORES

pdf(file="argleu_scores_by_kingdom_RGF_allGC.pdf", height=7, width=7)
plot(argleu_pc1_all_allGC, argleu_pc2_all_allGC, type="n", xlab="PC1", ylab="PC2")
title(main="ARG-LEU codon usage across kingdoms (all GC's)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(argleu_pc1_arch_allGC, argleu_pc2_arch_allGC, pch=15, col="red")
points(argleu_pc1_bact_allGC, argleu_pc2_bact_allGC, pch=16, col="purple")
points(argleu_pc1_euk_allGC, argleu_pc2_euk_allGC, pch=17,col="forestgreen")
#arrows(0, 0, X[,1], X[,2], len=0.1, col="black")
#text(1.4*X, rownames(X), col="black", xpd=T)
legend('topright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="argleu_scores_by_kingdom_RGF_below30.pdf", height=7, width=7)
plot(argleu_pc1_all_below30, argleu_pc2_all_below30, type="n", xlab="PC1", ylab="PC2")
title(main="ARG-LEU codon usage across kingdoms (GC < 30%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
#points(argleu_pc1_arch_below30, argleu_pc2_arch_below30, pch=15, col="red")
points(argleu_pc1_bact_below30, argleu_pc2_bact_below30, pch=16, col="purple")
points(argleu_pc1_euk_below30, argleu_pc2_euk_below30, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Bacteria','Eukarya'),col=c('purple','forestgreen'),pch=c(16,17), box.col="white")
dev.off()

pdf(file="argleu_scores_by_kingdom_RGF_30_40.pdf", height=7, width=7)
plot(argleu_pc1_all_30_40, argleu_pc2_all_30_40, type="n", xlab="PC1", ylab="PC2")
title(main="ARG-LEU codon usage across kingdoms (GC = 30-40%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(argleu_pc1_arch_30_40, argleu_pc2_arch_30_40, pch=15, col="red")
points(argleu_pc1_bact_30_40, argleu_pc2_bact_30_40, pch=16, col="purple")
points(argleu_pc1_euk_30_40, argleu_pc2_euk_30_40, pch=17,col="forestgreen")
legend('bottomright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="argleu_scores_by_kingdom_RGF_40_50.pdf", height=7, width=7)
plot(argleu_pc1_all_40_50, argleu_pc2_all_40_50, type="n", xlab="PC1", ylab="PC2")
title(main="ARG-LEU codon usage across kingdoms (GC = 40-50%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(argleu_pc1_arch_40_50, argleu_pc2_arch_40_50, pch=15, col="red")
points(argleu_pc1_bact_40_50, argleu_pc2_bact_40_50, pch=16, col="purple")
points(argleu_pc1_euk_40_50, argleu_pc2_euk_40_50, pch=17,col="forestgreen")
legend('bottomright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="argleu_scores_by_kingdom_RGF_50_60.pdf", height=7, width=7)
plot(argleu_pc1_all_50_60, argleu_pc2_all_50_60, type="n", xlab="PC1", ylab="PC2")
title(main="ARG-LEU codon usage across kingdoms (GC = 50-60%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(argleu_pc1_arch_50_60, argleu_pc2_arch_50_60, pch=15, col="red")
points(argleu_pc1_bact_50_60, argleu_pc2_bact_50_60, pch=16, col="purple")
points(argleu_pc1_euk_50_60, argleu_pc2_euk_50_60, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="argleu_scores_by_kingdom_RGF_60_70.pdf", height=7, width=7)
plot(argleu_pc1_all_60_70, argleu_pc2_all_60_70, type="n", xlab="PC1", ylab="PC2")
title(main="ARG-LEU codon usage across kingdoms (GC = 60-70%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(argleu_pc1_arch_60_70, argleu_pc2_arch_60_70, pch=15, col="red")
points(argleu_pc1_bact_60_70, argleu_pc2_bact_60_70, pch=16, col="purple")
points(argleu_pc1_euk_60_70, argleu_pc2_euk_60_70, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

####### REBUILD ANALYSIS FOR PC2-PC3 (only scores for ALL codons)

# PLOT PCA SCORES PC1-PC2

pdf(file="scores_by_kingdom_RGF_allGC_pc2pc3.pdf", height=7, width=7)
plot(pc2_all_allGC, pc3_all_allGC, type="n", xlab="PC2", ylab="PC3")
title(main="Codon usage across kingdoms (all GC's)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc2_arch_allGC, pc3_arch_allGC, pch=15, col="red")
points(pc2_bact_allGC, pc3_bact_allGC, pch=16, col="purple")
points(pc2_euk_allGC, pc3_euk_allGC, pch=17,col="forestgreen")
#arrows(0, 0, X[,1], X[,2], len=0.1, col="black")
#text(1.4*X, rownames(X), col="black", xpd=T)
legend('topright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="scores_by_kingdom_RGF_below30_pc2pc3.pdf", height=7, width=7)
plot(pc2_all_below30, pc3_all_below30, type="n", xlab="PC2", ylab="PC3")
title(main="Codon usage across kingdoms (GC < 30%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
#points(pc2_arch_below30, pc3_arch_below30, pch=15, col="red")
points(pc2_bact_below30, pc3_bact_below30, pch=16, col="purple")
points(pc2_euk_below30, pc3_euk_below30, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Bacteria','Eukarya'),col=c('purple','forestgreen'),pch=c(16,17), box.col="white")
dev.off()

pdf(file="scores_by_kingdom_RGF_30_40_pc2pc3.pdf", height=7, width=7)
plot(pc2_all_30_40, pc3_all_30_40, type="n", xlab="PC2", ylab="PC3")
title(main="Codon usage across kingdoms (GC = 30-40%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc2_arch_30_40, pc3_arch_30_40, pch=15, col="red")
points(pc2_bact_30_40, pc3_bact_30_40, pch=16, col="purple")
points(pc2_euk_30_40, pc3_euk_30_40, pch=17,col="forestgreen")
legend('bottomright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="scores_by_kingdom_RGF_40_50_pc2pc3.pdf", height=7, width=7)
plot(pc2_all_40_50, pc3_all_40_50, type="n", xlab="PC2", ylab="PC3")
title(main="Codon usage across kingdoms (GC = 40-50%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc2_arch_40_50, pc3_arch_40_50, pch=15, col="red")
points(pc2_bact_40_50, pc3_bact_40_50, pch=16, col="purple")
points(pc2_euk_40_50, pc3_euk_40_50, pch=17,col="forestgreen")
legend('bottomright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="scores_by_kingdom_RGF_50_60_pc2pc3.pdf", height=7, width=7)
plot(pc2_all_50_60, pc3_all_50_60, type="n", xlab="PC2", ylab="PC3")
title(main="Codon usage across kingdoms (GC = 50-60%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc2_arch_50_60, pc3_arch_50_60, pch=15, col="red")
points(pc2_bact_50_60, pc3_bact_50_60, pch=16, col="purple")
points(pc2_euk_50_60, pc3_euk_50_60, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()

pdf(file="scores_by_kingdom_RGF_60_70_pc2pc3.pdf", height=7, width=7)
plot(pc2_all_60_70, pc3_all_60_70, type="n", xlab="PC2", ylab="PC3")
title(main="Codon usage across kingdoms (GC = 60-70%)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc2_arch_60_70, pc3_arch_60_70, pch=15, col="red")
points(pc2_bact_60_70, pc3_bact_60_70, pch=16, col="purple")
points(pc2_euk_60_70, pc3_euk_60_70, pch=17,col="forestgreen")
legend('bottomleft', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()


