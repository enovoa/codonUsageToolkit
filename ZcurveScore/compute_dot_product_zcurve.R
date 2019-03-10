## COMPUTE DOT PRODUCT FROM PARAM FILE (including offset as last element) and ZVECTOR AND PLOT

arg <- commandArgs(trailingOnly = TRUE)

# 1. Read data
params<-read.table(arg[1]) 
testSetCoding<-read.table(arg[2])
print ("Read testSetCoding")
testSetNonCoding<-read.table(arg[3])
print ("Read testSetNonCoding")

fileTestSetCoding<-paste(arg[2],".Rdata",sep="")
fileTestSetNoncoding<-paste(arg[3],".Rdata",sep="")
save(testSetCoding, file=fileTestSetCoding)
save(testSetNonCoding, file=fileTestSetNoncoding)
print ("Saved RData files for testSetCoding and testSetNonCoding")

#example:


# 2. Compute dot product
a<-as.matrix(testSetCoding)
b<-as.matrix(testSetNonCoding)
p<-as.matrix(params[1:dim(params)[1]-1,])
res1<-a%*%p
res2<-b%*%p
zscores.testCoding<-res1+params[dim(params)[1],]
zscores.testNonCoding<-res2+params[dim(params)[1],]

# 3. Write zcurves
zcurveTestSetCodingFile<-paste("ZCURVE.",arg[2],sep="")
zcurveTestSetNonCodingFile<-paste("ZCURVE.",arg[3], sep="")

write.table(zscores.testCoding, file=zcurveTestSetCodingFile,quote=FALSE, sep="\t",  col.names=FALSE, row.names=FALSE)
write.table(zscores.testNonCoding, file=zcurveTestSetNonCodingFile,quote=FALSE, sep="\t",  col.names=FALSE, row.names=FALSE)

print ("Written tables with zcurve scores")

# 4. Plot
pdfName1<-paste(arg[2],".VS.",sep="")
pdfName2<-paste(pdfName1,arg[3],sep="")
pdfName<-paste(pdfName2,".pdf",sep="")
pdf(file=pdfName, height=7,width=7)
plot(density(zscores.testCoding[,1]), col="blue", xlab="ZcurveScore", ylim=c(0,1), xlim=c(-6,6), main="EUK VS BACT",lwd=4)
lines(density(zscores.testNonCoding[,1]), col="green",lwd=4)
legend('topleft',c(arg[2],arg[3]),lty=c(1,1),lwd=c(2.5,2.5),col=c("blue", "green"), box.col="white")
dev.off()                                  

print ("Written pdf plot of zscores") 

print ("Done! :-)")

