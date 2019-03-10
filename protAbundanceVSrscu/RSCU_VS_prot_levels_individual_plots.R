### PLOTTING PROT EXPRESSION VS CODON USAGE PER AMINO ACID

source("src/convert_RSCU_mat_into_codon_freq_df.R")
source("src/FUNCTIONS/RSCU_vs_prot_abundace_plot.R")

#### 1. Get data
# STRUCTURE OF THE DATASET:
# First row is prot name
# Codon usage RSCU for each 61 codons ("ter" removed)
# First row is AA_Codon column

raw_data<-read.table("dat/scerevisiae_RSCU_prot_expression_FINAL.txt", header=T)
dat<-raw_data[,2:339]
row.names(dat)<-raw_data$Codon

prot_levels<-read.table("dat/scerevisiae_ids_with_prot_levels.txt", header=T)
#First row is prot names
# Second row is prot values
row.names(prot_levels)<-"Prot_levels"

#### 2. Correlate both datasets

merged_dat<-rbind(dat,prot_levels)
final_dat<-as.data.frame(t(merged_dat))
codon_data<-final_dat[,1:61]
prot_data<-final_dat[,62]

# Check data well built
head(codon_data)[,1:5]
head(prot_data)

#### 3. Build plots for each amino acid

pdf(file="results/scerevisiae_plots_per_aa_RSCU_vs_prot_abundance_including_trendlines.pdf",height=7,width=7)

# a) 2 box
plot_RSCU_vs_prot_2box(final_dat$Lys_AAA,final_dat$Lys_AAG,final_dat$Prot_levels,final_dat,"Lys_AAA","Lys_AAG","Lys")
plot_RSCU_vs_prot_2box(final_dat$Asn_AAC,final_dat$Asn_AAU,final_dat$Prot_levels,final_dat,"Asn_AAC","Asn_AAU","Asn")
plot_RSCU_vs_prot_2box(final_dat$Gln_CAA,final_dat$Gln_CAG,final_dat$Prot_levels,final_dat,"Gln_CAA","Gln_CAG","Gln")
plot_RSCU_vs_prot_2box(final_dat$His_CAC,final_dat$His_CAU,final_dat$Prot_levels,final_dat,"His_CAC","His_CAU","His")
plot_RSCU_vs_prot_2box(final_dat$Glu_GAA,final_dat$Glu_GAG,final_dat$Prot_levels,final_dat,"Glu_GAA","Glu_GAG","Glu")
plot_RSCU_vs_prot_2box(final_dat$Asp_GAC,final_dat$Asp_GAU,final_dat$Prot_levels,final_dat,"Asp_GAC","Asp_GAU","Asp")
plot_RSCU_vs_prot_2box(final_dat$Tyr_UAU,final_dat$Tyr_UAC,final_dat$Prot_levels,final_dat,"Tyr_UAU","Tyr_UAC","Tyr")
plot_RSCU_vs_prot_2box(final_dat$Cys_UGU,final_dat$Cys_UGC,final_dat$Prot_levels,final_dat,"Cys_UGU","Cys_UGC","Cys")
plot_RSCU_vs_prot_2box(final_dat$Phe_UUU,final_dat$Phe_UUC,final_dat$Prot_levels,final_dat,"Phe_UUU","Phe_UUC","Phe")

# b) 4 box
plot_RSCU_vs_prot_4box(final_dat$Ala_GCA,final_dat$Ala_GCC, final_dat$Ala_GCG, final_dat$Ala_GCU,final_dat$Prot_levels,final_dat,"Ala_GCA","Ala_GCC","Ala_GCG","Ala_GCU","Ala")
plot_RSCU_vs_prot_4box(final_dat$Thr_ACA,final_dat$Thr_ACC, final_dat$Thr_ACG, final_dat$Thr_ACU,final_dat$Prot_levels,final_dat,"Thr_ACA","Thr_ACC","Thr_ACG","Thr_ACU","Thr")
plot_RSCU_vs_prot_4box(final_dat$Gly_GGA,final_dat$Gly_GGC, final_dat$Gly_GGG, final_dat$Gly_GGU,final_dat$Prot_levels,final_dat,"Gly_GGA","Gly_GGC","Gly_GGG","Gly_GGU","Gly")
plot_RSCU_vs_prot_4box(final_dat$Val_GUA,final_dat$Val_GUC, final_dat$Val_GUG, final_dat$Val_GUU,final_dat$Prot_levels,final_dat,"Val_GUA","Val_GUC","Val_GUG","Val_GUU","Val")
plot_RSCU_vs_prot_4box(final_dat$Pro_CCA,final_dat$Pro_CCC, final_dat$Pro_CCG, final_dat$Pro_CCU,final_dat$Prot_levels,final_dat,"Pro_CCA","Pro_CCC","Pro_CCG","Pro_CCU","Pro")

# c) 3 box
plot_RSCU_vs_prot_3box(final_dat$Ile_AUA,final_dat$Ile_AUC, final_dat$Ile_AUU,final_dat$Prot_levels,final_dat,"Ile_AUA","Ile_AUC","Ile_AUU","Ile")

# d) 6 box
plot_RSCU_vs_prot_6box(final_dat$Ser_UCA,final_dat$Ser_UCC, final_dat$Ser_UCG, final_dat$Ser_UCU, final_dat$Ser_AGC, final_dat$Ser_AGU, final_dat$Prot_levels,final_dat,"Ser_UCA","Ser_UCC","Ser_UCG","Ser_UCU","Ser_AGC","Ser_AGU","Ser")
plot_RSCU_vs_prot_6box(final_dat$Leu_CUA,final_dat$Leu_CUC, final_dat$Leu_CUG, final_dat$Leu_CUU, final_dat$Leu_UUA, final_dat$Leu_UUG, final_dat$Prot_levels,final_dat,"Leu_CUA","Leu_CUC","Leu_CUG","Leu_CUU","Leu_UUA","Leu_UUG","Leu")
plot_RSCU_vs_prot_6box(final_dat$Arg_CGA,final_dat$Arg_CGC, final_dat$Arg_CGG, final_dat$Arg_CGU, final_dat$Arg_AGA, final_dat$Arg_AGG, final_dat$Prot_levels,final_dat,"Arg_CGA","Arg_CGC","Arg_CGG","Arg_CGU","Arg_AGA","Arg_AGG","Arg")


dev.off()



