# 1. READ RSCU DATA PER-SPECIES

dat_allGC_with_species<-read.table("dat/all_allGC.txt.transposed", header=T) #dat_allGC, ordered RSCU values, 120 first species are ARCH, next 1341 are BACT, next 164 are EUK
row.names(dat_allGC_with_species)<-paste(dat_allGC_with_species$Species,row.names(dat_allGC_with_species),sep=".")
dat_allGC<-dat_allGC_with_species[,2:65]

# 2. BUILD PCA
model_allGC<-prcomp(dat_allGC)
summary(model_allGC)

# 3. GET SCORES MATRIX PCA
scores_allGC<-model_allGC$x
scores_data_allGC<-as.matrix(scores_allGC[,1:3])
pc1_all_allGC<-as.vector(scores_data_allGC[,1])
pc2_all_allGC<-as.vector(scores_data_allGC[,2])
pc3_all_allGC<-as.vector(scores_data_allGC[,3])

# 4. SUBDIVIDE SCORES PER KINGDOM
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

# 5. GET LOADINGS MATRIX PCA
loadings_allGC<-model_allGC$rotation
loadings_data_allGC<-as.matrix(loadings_allGC[,1:3])
loadings_pc1_all_allGC<-as.vector(loadings_data_allGC[,1])
loadings_pc2_all_allGC<-as.vector(loadings_data_allGC[,2])
loadings_pc3_all_allGC<-as.vector(loadings_data_allGC[,3])

# 6. PLOT PCA SCORES PC1-PC2

#pdf(file="scores_by_kingdom_RGF_allGC.pdf", height=7, width=7)
plot(pc1_all_allGC, pc2_all_allGC, type="n", xlab="PC1", ylab="PC2")
title(main="Codon usage across kingdoms (all GC's)", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc1_arch_allGC, pc2_arch_allGC, pch=15, col="red")
points(pc1_bact_allGC, pc2_bact_allGC, pch=16, col="purple")
points(pc1_euk_allGC, pc2_euk_allGC, pch=17,col="forestgreen")
legend('topright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
#dev.off()

# 7. PLOT PCA LOADINGS PC1-PC2

# Separate labels in order to color according to its last codon --> use data_allGC, which can be used for all plots
labels_A<-cbind(rownames(loadings_data_allGC)[1],rownames(loadings_data_allGC)[5],rownames(loadings_data_allGC)[9],rownames(loadings_data_allGC)[13],rownames(loadings_data_allGC)[17],rownames(loadings_data_allGC)[21],rownames(loadings_data_allGC)[25],rownames(loadings_data_allGC)[29],rownames(loadings_data_allGC)[33],rownames(loadings_data_allGC)[37],rownames(loadings_data_allGC)[41],rownames(loadings_data_allGC)[45],rownames(loadings_data_allGC)[49],rownames(loadings_data_allGC)[53],rownames(loadings_data_allGC)[57],rownames(loadings_data_allGC)[61])
labels_C<-cbind(rownames(loadings_data_allGC)[2],rownames(loadings_data_allGC)[6],rownames(loadings_data_allGC)[10],rownames(loadings_data_allGC)[14],rownames(loadings_data_allGC)[18],rownames(loadings_data_allGC)[22],rownames(loadings_data_allGC)[26],rownames(loadings_data_allGC)[30],rownames(loadings_data_allGC)[34],rownames(loadings_data_allGC)[38],rownames(loadings_data_allGC)[42],rownames(loadings_data_allGC)[46],rownames(loadings_data_allGC)[50],rownames(loadings_data_allGC)[54],rownames(loadings_data_allGC)[58],rownames(loadings_data_allGC)[62])
labels_G<-cbind(rownames(loadings_data_allGC)[3],rownames(loadings_data_allGC)[7],rownames(loadings_data_allGC)[11],rownames(loadings_data_allGC)[15],rownames(loadings_data_allGC)[19],rownames(loadings_data_allGC)[23],rownames(loadings_data_allGC)[27],rownames(loadings_data_allGC)[31],rownames(loadings_data_allGC)[35],rownames(loadings_data_allGC)[39],rownames(loadings_data_allGC)[43],rownames(loadings_data_allGC)[47],rownames(loadings_data_allGC)[51],rownames(loadings_data_allGC)[55],rownames(loadings_data_allGC)[59],rownames(loadings_data_allGC)[63])
labels_U<-cbind(rownames(loadings_data_allGC)[4],rownames(loadings_data_allGC)[8],rownames(loadings_data_allGC)[12],rownames(loadings_data_allGC)[16],rownames(loadings_data_allGC)[20],rownames(loadings_data_allGC)[24],rownames(loadings_data_allGC)[28],rownames(loadings_data_allGC)[32],rownames(loadings_data_allGC)[36],rownames(loadings_data_allGC)[40],rownames(loadings_data_allGC)[44],rownames(loadings_data_allGC)[48],rownames(loadings_data_allGC)[52],rownames(loadings_data_allGC)[56],rownames(loadings_data_allGC)[60],rownames(loadings_data_allGC)[64])


#pdf(file="loadings_by_kingdom_RGF_allGC.pdf", height=7, width=7)
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
#dev.off()


##########################

# ANALYSIS OF ONLY ARG ###

##########################

arg_dat_allGC<-cbind(dat_allGC$CGU, dat_allGC$CGC, dat_allGC$CGA, dat_allGC$CGG, dat_allGC$AGA, dat_allGC$AGG)
colnames(arg_dat_allGC)=c("CGU","CGC","CGA","CGG","AGA","AGG")
arg_model_allGC<-prcomp(arg_dat_allGC)
arg_scores_allGC<-arg_model_allGC$x
arg_scores_data_allGC<-as.matrix(arg_scores_allGC[,1:2])
arg_pc1_all_allGC<-as.vector(arg_scores_data_allGC[,1])
arg_pc2_all_allGC<-as.vector(arg_scores_data_allGC[,2])

#GET SCORES FOR ARCH,EUK AND BACT

arg_arch_allGC<-arg_scores_data_allGC[1:120,] 		# 120
arg_bact_allGC<-arg_scores_data_allGC[121:1461,]  	# 1341
arg_euk_allGC<-arg_scores_data_allGC[1462:1625,]	# 164 --> total = 1625
arg_pc1_arch_allGC<-as.vector(arg_arch_allGC[,1])
arg_pc2_arch_allGC<-as.vector(arg_arch_allGC[,2])
arg_pc1_euk_allGC<-as.vector(arg_euk_allGC[,1])
arg_pc2_euk_allGC<-as.vector(arg_euk_allGC[,2])
arg_pc1_bact_allGC<-as.vector(arg_bact_allGC[,1])
arg_pc2_bact_allGC<-as.vector(arg_bact_allGC[,2])


# PLOT PCA SCORES

#pdf(file="arg_scores_by_kingdom_RGF_allGC.pdf", height=7, width=7)
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
#dev.off()


#################################

# ANALYSIS OF  ALL EXCEPT ARG ###

#################################

arg_codons<-c("CGU","CGC","CGA","CGG","AGA","AGG")
noarg_dat_allGC<-dat_allGC[,!(colnames(dat_allGC) %in% arg_codons), drop = FALSE]

noarg_model_allGC<-prcomp(noarg_dat_allGC)
noarg_scores_allGC<-noarg_model_allGC$x
noarg_scores_data_allGC<-as.matrix(noarg_scores_allGC[,1:2])
noarg_pc1_all_allGC<-as.vector(noarg_scores_data_allGC[,1])
noarg_pc2_all_allGC<-as.vector(noarg_scores_data_allGC[,2])

#GET SCORES FOR ARCH,EUK AND BACT

noarg_arch_allGC<-noarg_scores_data_allGC[1:120,] 		# 120
noarg_bact_allGC<-noarg_scores_data_allGC[121:1461,]  	# 1341
noarg_euk_allGC<-noarg_scores_data_allGC[1462:1625,]	# 164 --> total = 1625
noarg_pc1_arch_allGC<-as.vector(noarg_arch_allGC[,1])
noarg_pc2_arch_allGC<-as.vector(noarg_arch_allGC[,2])
noarg_pc1_bact_allGC<-as.vector(noarg_bact_allGC[,1])
noarg_pc2_bact_allGC<-as.vector(noarg_bact_allGC[,2])
noarg_pc1_euk_allGC<-as.vector(noarg_euk_allGC[,1])
noarg_pc2_euk_allGC<-as.vector(noarg_euk_allGC[,2])

# PLOT PCA SCORES

#pdf(file="noarg_scores_by_kingdom_RGF_allGC.pdf", height=7, width=7)
plot(noarg_pc1_all_allGC, noarg_pc2_all_allGC, type="n", xlab="PC1", ylab="PC2")
title(main="noarg codon usage across kingdoms", col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(noarg_pc1_arch_allGC, noarg_pc2_arch_allGC, pch=15, col="red")
points(noarg_pc1_bact_allGC, noarg_pc2_bact_allGC, pch=16, col="purple")
points(noarg_pc1_euk_allGC, noarg_pc2_euk_allGC, pch=17,col="forestgreen")
legend('topleft', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
#dev.off()



