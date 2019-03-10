# 1. READ RSCU DATA PER-SPECIES

dat_allGC<-read.table("dat/per_species_codonusage_RSCU_ALL_emblcds.txt", header=T)# ordered RSCU values, 120 first species are ARCH, next 1341 are BACT, next 164 are EUK

# 2. BUILD PCA
model_allGC<-prcomp(dat_allGC)
summary(model_allGC) # PC1 = 71.6%; PC2 = 7.6%

# 3. GET SCORES MATRIX PCA
scores_allGC<-model_allGC$x
scores_data_allGC<-as.matrix(scores_allGC[,1:4])
pc1_all_allGC<-as.vector(scores_data_allGC[,1])
pc2_all_allGC<-as.vector(scores_data_allGC[,2])
pc3_all_allGC<-as.vector(scores_data_allGC[,3])
pc4_all_allGC<-as.vector(scores_data_allGC[,4])

# 4. SUBDIVIDE SCORES PER KINGDOM
arch_allGC<-scores_data_allGC[1:120,] 		# 120
bact_allGC<-scores_data_allGC[121:1461,]  	# 1341
euk_allGC<-scores_data_allGC[1462:1625,]	# 164 --> total = 1625
pc1_arch_allGC<-as.vector(arch_allGC[,1])
pc2_arch_allGC<-as.vector(arch_allGC[,2])
pc3_arch_allGC<-as.vector(arch_allGC[,3])
pc4_arch_allGC<-as.vector(arch_allGC[,4])
pc1_euk_allGC<-as.vector(euk_allGC[,1])
pc2_euk_allGC<-as.vector(euk_allGC[,2])
pc3_euk_allGC<-as.vector(euk_allGC[,3])
pc4_euk_allGC<-as.vector(euk_allGC[,4])
pc1_bact_allGC<-as.vector(bact_allGC[,1])
pc2_bact_allGC<-as.vector(bact_allGC[,2])
pc3_bact_allGC<-as.vector(bact_allGC[,3])
pc4_bact_allGC<-as.vector(bact_allGC[,4])

# 5. GET LOADINGS MATRIX PCA
loadings_allGC<-model_allGC$rotation
loadings_data_allGC<-as.matrix(loadings_allGC[,1:4])
loadings_pc1<-as.vector(loadings_data_allGC[,1])
loadings_pc2<-as.vector(loadings_data_allGC[,2])
loadings_pc3<-as.vector(loadings_data_allGC[,3])
loadings_pc4<-as.vector(loadings_data_allGC[,4])

# 6. PLOT PCA SCORES PC1-PC2

setwd("./results")

plot_pca_scores<-function(pc1_all_allGC,pc2_all_allGC,pc1_arch_allGC,pc2_arch_allGC,pc1_bact_allGC,pc2_bact_allGC,pc1_euk_allGC,pc2_euk_allGC,xlab,ylab,main) {
pdf(file=paste(main,".pdf",sep=""), height=7, width=7)
plot(pc1_all_allGC, pc2_all_allGC, type="n", xlab=xlab, ylab=ylab)
title(main=main, col.main="black", font.main=4)
abline(v=0, lty=3)
abline(h=0, lty=3)
points(pc1_arch_allGC, pc2_arch_allGC, pch=15, col="red")
points(pc1_bact_allGC, pc2_bact_allGC, pch=16, col="purple")
points(pc1_euk_allGC, pc2_euk_allGC, pch=17,col="forestgreen")
#legend('topright', cex=1.2,ncol=1,c('Archaea','Bacteria','Eukarya'),col=c('red','purple','forestgreen'),pch=c(15,16,17), box.col="white")
dev.off()
}

plot_pca_scores(pc1_all_allGC,pc2_all_allGC,pc1_arch_allGC,pc2_arch_allGC,pc1_bact_allGC,pc2_bact_allGC,pc1_euk_allGC,pc2_euk_allGC,"PC1","PC2","per-species codon usage across kingdoms-PC1-PC2")
plot_pca_scores(pc1_all_allGC,pc3_all_allGC,pc1_arch_allGC,pc3_arch_allGC,pc1_bact_allGC,pc3_bact_allGC,pc1_euk_allGC,pc3_euk_allGC,"PC1","PC3","per-species codon usage across kingdoms-PC1-PC3")
plot_pca_scores(pc1_all_allGC,pc4_all_allGC,pc1_arch_allGC,pc4_arch_allGC,pc1_bact_allGC,pc4_bact_allGC,pc1_euk_allGC,pc4_euk_allGC,"PC1","PC4","per-species codon usage across kingdoms-PC1-PC4")
plot_pca_scores(pc2_all_allGC,pc4_all_allGC,pc2_arch_allGC,pc4_arch_allGC,pc2_bact_allGC,pc4_bact_allGC,pc2_euk_allGC,pc4_euk_allGC,"PC2","PC4","per-species codon usage across kingdoms-PC2-PC4")
plot_pca_scores(pc2_all_allGC,pc3_all_allGC,pc2_arch_allGC,pc3_arch_allGC,pc2_bact_allGC,pc3_bact_allGC,pc2_euk_allGC,pc3_euk_allGC,"PC2","PC3","per-species codon usage across kingdoms-PC2-PC3")


# 7. PLOT PCA LOADINGS PC1-PC2

# Separate labels in order to color according to its last codon --> use data_allGC, which can be used for all plots
labels_A<-rownames(loadings_data_allGC)[c(1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61)]
labels_C<-rownames(loadings_data_allGC)[c(2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62)]
labels_G<-rownames(loadings_data_allGC)[c(3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63)]
labels_U<-rownames(loadings_data_allGC)[c(4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64)]


plot_loadings<-function(loadings_pc1,loadings_pc2,xlab,ylab, main) {
    # loadings pc1
    loadings_pc1_A<-loadings_pc1[c(1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61)]
    loadings_pc1_C<-loadings_pc1[c(2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62)]
    loadings_pc1_G<-loadings_pc1[c(3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63)]
    loadings_pc1_U<-loadings_pc1[c(4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64)]

    # loadings pc2
    loadings_pc2_A<-loadings_pc2[c(1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61)]
    loadings_pc2_C<-loadings_pc2[c(2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62)]
    loadings_pc2_G<-loadings_pc2[c(3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63)]
    loadings_pc2_U<-loadings_pc2[c(4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64)]

    pdf(file=paste(main,".pdf",sep=""), height=7, width=7)
    plot(loadings_pc1,loadings_pc2, type="n", xlab=xlab, ylab=ylab)
    title(main=main, col.main="black", font.main=4)
    abline(v=0, lty=3)
    abline(h=0, lty=3)
    text(loadings_pc1_A,loadings_pc2_A, labels_A, col="blue", xpd=T)
    text(loadings_pc1_C,loadings_pc2_C, labels_C, col="red", xpd=T)
    text(loadings_pc1_G,loadings_pc2_G, labels_G, col="orange", xpd=T)
    text(loadings_pc1_U,loadings_pc2_U, labels_U, col="purple", xpd=T)
    legend('bottomright', cex=1,ncol=1,c('A-ended','U-ended','G-ended','C-ended'),col=c('blue','purple','orange','red'),pch=c(15,15,15,15), box.col="white")
    dev.off()
}

plot_loadings(loadings_pc1,loadings_pc2,"PC1","PC2", "loadings_per_species-PC1-PC2")
plot_loadings(loadings_pc2,loadings_pc3,"PC2","PC3", "loadings_per_species-PC2-PC3")
plot_loadings(loadings_pc1,loadings_pc3,"PC1","PC3", "loadings_per_species-PC1-PC3")
plot_loadings(loadings_pc2,loadings_pc4,"PC2","PC4", "loadings_per_species-PC2-PC4")
plot_loadings(loadings_pc1,loadings_pc4,"PC1","PC4", "loadings_per_species-PC1-PC4")


##########################

# ANALYSIS OF ONLY ARG ###

##########################

arg_dat_allGC<-cbind(dat_allGC$CGU, dat_allGC$CGC, dat_allGC$CGA, dat_allGC$CGG, dat_allGC$AGA, dat_allGC$AGG)
colnames(arg_dat_allGC)=c("CGU","CGC","CGA","CGG","AGA","AGG")
arg_model_allGC<-prcomp(arg_dat_allGC)
summary(arg_model_allGC)
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

pdf(file="arg_scores_by_kingdom.pdf", height=7, width=7)
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

# BOXPLOTS
pdf(file="arg_boxplot.pdf", height=7, width=7)
arg_dat_allGC<-arg_dat_allGC[,order(colnames(arg_dat_allGC))]
arg_dat_arch<-arg_dat_allGC[1:120,] 		# 120
arg_dat_bact<-arg_dat_allGC[121:1461,]  	# 1341
arg_dat_euk<-arg_dat_allGC[1462:1625,]	# 164 --> total = 1625
par(mfrow=c(3,1))
boxplot(arg_dat_arch, col="red",cex.axis=2)
abline(h=1,lty=2)
boxplot(arg_dat_bact, col="purple",cex.axis=2)
abline(h=1,lty=2)
boxplot(arg_dat_euk, col="forestgreen",cex.axis=2)
abline(h=1,lty=2)
dev.off()

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



