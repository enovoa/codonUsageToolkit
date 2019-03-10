### AVERAGE CODON USAGE PER SPECIES - including DOMAIN and GC CONTENT info  ###


## 1. READ DATA (by kingdom and GC content) ##
directory
arch_30_40<-read.table("dat/arch_30-40.pasted.sorted.values.transposed")
arch_40_50<-read.table("dat/arch_40-50.pasted.sorted.values.transposed")
arch_50_60<-read.table("dat/arch_50-60.pasted.sorted.values.transposed")
arch_60_70<-read.table("dat/arch_60-70.pasted.sorted.values.transposed")
arch_below30<-read.table("dat/arch_below30.pasted.sorted.values.transposed")
arch_over70<-read.table("dat/arch_over70.pasted.sorted.values.transposed")

bact_30_40<-read.table("dat/bact_30-40.pasted.sorted.values.transposed")
bact_40_50<-read.table("dat/bact_40-50.pasted.sorted.values.transposed")
bact_50_60<-read.table("dat/bact_50-60.pasted.sorted.values.transposed")
bact_60_70<-read.table("dat/bact_60-70.pasted.sorted.values.transposed")
bact_below30<-read.table("dat/bact_below30.pasted.sorted.values.transposed")
bact_over70<-read.table("dat/bact_over70.pasted.sorted.values.transposed")

euk_30_40<-read.table("dat/euk_30-40.pasted.sorted.values.transposed")
euk_40_50<-read.table("dat/euk_40-50.pasted.sorted.values.transposed")
euk_50_60<-read.table("dat/euk_50-60.pasted.sorted.values.transposed")
euk_60_70<-read.table("dat/euk_60-70.pasted.sorted.values.transposed")
euk_below30<-read.table("dat/euk_below30.pasted.sorted.values.transposed")
euk_over70<-read.table("dat/euk_over70.pasted.sorted.values.transposed")

# Define GC vectors per kingdom
gc_arch.group<-c(rep("yellow",dim(arch_below30)[1]), rep("orange",dim(arch_30_40)[1]),rep("red",dim(arch_40_50)[1]),rep("darkred",dim(arch_50_60)[1]),rep("brown",dim(arch_60_70)[1]))
gc_bact.group<-c(rep("yellow",dim(bact_below30)[1]), rep("orange",dim(bact_30_40)[1]),rep("red",dim(bact_40_50)[1]),rep("darkred",dim(bact_50_60)[1]),rep("brown",dim(bact_60_70)[1]),rep("black",dim(bact_over70)[1]))
gc_euk.group<-c(rep("yellow",dim(euk_below30)[1]), rep("orange",dim(euk_30_40)[1]),rep("red",dim(euk_40_50)[1]),rep("darkred",dim(euk_50_60)[1]),rep("brown",dim(euk_60_70)[1]),rep("black",dim(euk_over70)[1]))
all_gc.group<-c(gc_arch.group,gc_bact.group,gc_euk.group)


# Merged dataset
all_arch<-rbind(arch_below30,arch_30_40,arch_40_50,arch_50_60,arch_60_70) # 120 sps
all_bact<-rbind(bact_below30,bact_30_40,bact_40_50,bact_50_60,bact_60_70,bact_over70) # 1341 sps
all_euk<-rbind(euk_below30,euk_30_40,euk_40_50,euk_50_60,euk_60_70,euk_over70) #164 sps
all_allkingdoms<-rbind(all_arch, all_bact, all_euk)

# Define kingdom vector
kingdom.group<-c(rep("blue",dim(all_arch)[1]),rep("darkgoldenrod1",dim(all_bact)[1]),rep("green4 ",dim(all_euk)[1]))

# Read codons
aa_codon<-read.table("dat/aa_and_codon_in_datasets.txt", header=T)

## 2. HEATMAP  ####

library(gplots)

# Distance and clustering functions
mydist=function(c) {dist(c,method="manhattan")}
myclust=function(c) {hclust(c,method="complete")}

# Define breaks
color.breaks <- c(seq(0, 0.79, length=100), 			# for color 1
				seq(0.8, 1.2, length=100),				# for color 2
				seq(1.21, 6, length=100))				# for color 3

# Define color palette
mycol <- colorpanel(n=299,low="green",mid="black",high="red")
#my_palette <- colorRampPalette(c("green","black","red"))(n=299)

# Add codon labels
colnames(all_allkingdoms)<-aa_codon$AA_Codon

# Remove ter, trp, met
d1<-all_allkingdoms[,-59]
d2<-d1[,-57]
d3<-d2[,-51]
d4<-d3[,-49]
all_allkingdoms.without_ter_trp_met<-d4[,-15]

# Plot
#png("heatmap_1625sps_average_codon_usage_GCcontent.png",
	width = 5 * 300,
	height = 5*300,
	res = 300,
	pointsize = 8 )

# Heatmap
heatmap.2(as.matrix(t(all_allkingdoms.without_ter_trp_met)),
	#cellnote = as.matrix(dat2), 	# cell labels with numbers
	density.info="none",			# turns off density plot inside color legend
	trace="none",					# turns off trace lines inside the heatmap
	key=TRUE, 						# show key			
	symkey=FALSE,					# key does not have to be symmetrical
	scale="none", 					# do not scale rows or columns
	cexCol=0.7, 					# font for column labels
	cexRow=0.7, 					# font for row labels
	col=mycol,						# color scheme
    breaks=color.breaks, 			# breaks for colors
	#labRow = "",					# row label
	labCol = "",					# column label
	#dendrogram = "row",				# only draw row dendogram
	#Colv = "NA",					# turn off column clustering
	#Rowv = "NA",
	ColSideColors=kingdom.group,		# side column coloring (a given characteristic chosen)
	#ColSideColors=all_gc.group,
	hclustfun=myclust, 				# clustering method
	distfun=mydist) 				# distance method

#dev.off()


