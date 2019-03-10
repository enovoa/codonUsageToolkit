# Create plot for RSCU for each amino acid VS protein abundance including loess trend line (local polynomial regression)

plot_RSCU_vs_prot_2box<- function(x1,x2,y,z,codon1,codon2,aa) {
	# x1,x2 is the data$codon
	# y is the prot abundance
	# z is the dataset that includes x1,x2 and y
	# codon1, codon2 are the codon names (for legend)
	# aa is the aa name (for main title of plot)
	
	plot(log(y), x1, col=rgb(0,0,1,1/4), ylim=c(0,2.3), pch=16, ylab="", xlab="", main=aa, cex.axis=1.3, cex.main=2)
	mtext("log(protein abundance)", side=1, line=3, cex=1.2)
	mtext("RSCU", side=2, line=3, cex=1.2)

	# fit a loess line
	loess_fit1 <- loess(x1~log(y), data=z)
	dat1<-as.data.frame(cbind(log(y),predict(loess_fit1)))
	colnames(dat1)<-c("prot","loess")
	dat1<-dat1[order(dat1$prot),]
	lines(dat1$prot, dat1$loess, col="blue", lwd=2)

	points(log(y), x2, col=rgb(1,0,0,1/4), pch=16)
	loess_fit2 <- loess(x2~log(y), data=z)
	dat2<-as.data.frame(cbind(log(y),predict(loess_fit2)))
	colnames(dat2)<-c("prot","loess")
	dat2<-dat2[order(dat2$prot),]
	lines(dat2$prot, dat2$loess, col="red", lwd=2)

	legend('topleft',c(codon1,codon2),lty=c(1,1),lwd=c(2.5,2.5),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),bty = "n", cex=1.3)
}

plot_RSCU_vs_prot_4box<- function(x1,x2,x3,x4,y,z,codon1,codon2,codon3,codon4,aa) {
	# x1,x2 is the data$codon
	# y is the prot abundance
	# z is the dataset that includes x1,x2 and y
	# codon1, codon2 are the codon names (for legend)
	# aa is the aa name (for main title of plot)
	
	plot(log(y), x1, col=rgb(0,0,1,1/4), ylim=c(0,4.6), pch=16, ylab="", xlab="", main=aa, cex.axis=1.3, cex.main=2)
	mtext("log(protein abundance)", side=1, line=3, cex=1.2)
	mtext("RSCU", side=2, line=3, cex=1.2)
	# fit a loess line
	loess_fit1 <- loess(x1~log(y), data=z)
	dat1<-as.data.frame(cbind(log(y),predict(loess_fit1)))
	colnames(dat1)<-c("prot","loess")
	dat1<-dat1[order(dat1$prot),]
	lines(dat1$prot, dat1$loess, col="blue", lwd=2)
	
	points(log(y), x2, col=rgb(1,0,0,1/4), pch=16)
	loess_fit2 <- loess(x2~log(y), data=z)
	dat2<-as.data.frame(cbind(log(y),predict(loess_fit2)))
	colnames(dat2)<-c("prot","loess")
	dat2<-dat2[order(dat2$prot),]
	lines(dat2$prot, dat2$loess, col="red", lwd=2)
	
	points(log(y), x3, col=rgb(0,1,1,1/4), pch=16)
	loess_fit3 <- loess(x3~log(y), data=z)
	dat3<-as.data.frame(cbind(log(y),predict(loess_fit3)))
	colnames(dat3)<-c("prot","loess")
	dat3<-dat3[order(dat3$prot),]
	lines(dat3$prot, dat3$loess, col="cyan", lwd=2)
	
	points(log(y), x4, col=rgb(0,1,0,1/4), pch=16)
	loess_fit4 <- loess(x4~log(y), data=z)
	dat4<-as.data.frame(cbind(log(y),predict(loess_fit4)))
	colnames(dat4)<-c("prot","loess")
	dat4<-dat4[order(dat4$prot),]
	lines(dat4$prot, dat4$loess, col="green", lwd=2)
	legend('topleft',c(codon1,codon2,codon3,codon4),lty=c(1,1),lwd=c(2.5,2.5),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4),rgb(0,1,1,1/4),rgb(0,1,0,1/4)), ncol=2,bty = "n", cex=1.3)
}

plot_RSCU_vs_prot_3box<- function(x1,x2,x3,y,z,codon1,codon2,codon3,aa) {
	# x1,x2 is the data$codon
	# y is the prot abundance
	# z is the dataset that includes x1,x2 and y
	# codon1, codon2 are the codon names (for legend)
	# aa is the aa name (for main title of plot)
	
	plot(log(y), x1, col=rgb(0,0,1,1/4), ylim=c(0,3.5), pch=16, ylab="", xlab="", main=aa, cex.axis=1.3, cex.main=2)
	mtext("log(protein abundance)", side=1, line=3, cex=1.2)
	mtext("RSCU", side=2, line=3, cex=1.2)
	# fit a loess line
	loess_fit1 <- loess(x1~log(y), data=z)
	dat1<-as.data.frame(cbind(log(y),predict(loess_fit1)))
	colnames(dat1)<-c("prot","loess")
	dat1<-dat1[order(dat1$prot),]
	lines(dat1$prot, dat1$loess, col="blue", lwd=2)
	
	points(log(y), x2, col=rgb(1,0,0,1/4), pch=16)
	loess_fit2 <- loess(x2~log(y), data=z)
	dat2<-as.data.frame(cbind(log(y),predict(loess_fit2)))
	colnames(dat2)<-c("prot","loess")
	dat2<-dat2[order(dat2$prot),]
	lines(dat2$prot, dat2$loess, col="red", lwd=2)
	
	points(log(y), x3, col=rgb(0,1,1,1/4), pch=16)
	loess_fit3 <- loess(x3~log(y), data=z)
	dat3<-as.data.frame(cbind(log(y),predict(loess_fit3)))
	colnames(dat3)<-c("prot","loess")
	dat3<-dat3[order(dat3$prot),]
	lines(dat3$prot, dat3$loess, col="cyan", lwd=2)
	
	legend('topleft',c(codon1,codon2,codon3),lty=c(1,1),lwd=c(2.5,2.5),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4),rgb(0,1,1,1/4)), ncol=2,bty = "n", cex=1.3)
}


plot_RSCU_vs_prot_6box<- function(x1,x2,x3,x4,x5,x6,y,z,codon1,codon2,codon3,codon4,codon5,codon6,aa) {
	# x1,x2 is the data$codon
	# y is the prot abundance
	# z is the dataset that includes x1,x2 and y
	# codon1, codon2 are the codon names (for legend)
	# aa is the aa name (for main title of plot)
	
	plot(log(y), x1, col=rgb(0,0,1,1/4), ylim=c(0,6.8), pch=16, ylab="", xlab="", main=aa, cex.axis=1.3, cex.main=2)
	mtext("log(protein abundance)", side=1, line=3, cex=1.2)
	mtext("RSCU", side=2, line=3, cex=1.2)
	# fit a loess line
	loess_fit1 <- loess(x1~log(y), data=z)
	dat1<-as.data.frame(cbind(log(y),predict(loess_fit1)))
	colnames(dat1)<-c("prot","loess")
	dat1<-dat1[order(dat1$prot),]
	lines(dat1$prot, dat1$loess, col="blue", lwd=2)
	
	points(log(y), x2, col=rgb(1,0,0,1/4), pch=16)
	loess_fit2 <- loess(x2~log(y), data=z)
	dat2<-as.data.frame(cbind(log(y),predict(loess_fit2)))
	colnames(dat2)<-c("prot","loess")
	dat2<-dat2[order(dat2$prot),]
	lines(dat2$prot, dat2$loess, col="red", lwd=2)
	
	points(log(y), x3, col=rgb(0,1,1,1/4), pch=16)
	loess_fit3 <- loess(x3~log(y), data=z)
	dat3<-as.data.frame(cbind(log(y),predict(loess_fit3)))
	colnames(dat3)<-c("prot","loess")
	dat3<-dat3[order(dat3$prot),]
	lines(dat3$prot, dat3$loess, col="cyan", lwd=2)
	
	points(log(y), x4, col=rgb(0,1,0,1/4), pch=16)
	loess_fit4 <- loess(x4~log(y), data=z)
	dat4<-as.data.frame(cbind(log(y),predict(loess_fit4)))
	colnames(dat4)<-c("prot","loess")
	dat4<-dat4[order(dat4$prot),]
	lines(dat4$prot, dat4$loess, col="green", lwd=2)
	
	points(log(y), x5, col=rgb(1,0,1,1/4), pch=16)
	loess_fit5 <- loess(x5~log(y), data=z)
	dat5<-as.data.frame(cbind(log(y),predict(loess_fit5)))
	colnames(dat5)<-c("prot","loess")
	dat5<-dat5[order(dat5$prot),]
	lines(dat5$prot, dat5$loess, col="magenta", lwd=2)
	
	points(log(y), x6, col=rgb(1,0.5,0,1/4), pch=16)
	loess_fit6 <- loess(x6~log(y), data=z)
	dat6<-as.data.frame(cbind(log(y),predict(loess_fit6)))
	colnames(dat6)<-c("prot","loess")
	dat6<-dat6[order(dat6$prot),]
	lines(dat6$prot, dat6$loess, col="orange", lwd=2)

	legend('topleft',c(codon1,codon2,codon3,codon4,codon5,codon6),lty=c(1,1),lwd=c(2.5,2.5),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4),rgb(0,1,1,1/4),rgb(0,1,0,1/4),rgb(1,0,1,1/4),rgb(1,0.5,0,1/4)), ncol=3,bty = "n", cex=1.3)
}





