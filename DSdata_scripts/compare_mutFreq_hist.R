pdf(file="compare_mutFreq_hist.pdf", width=8, height=8)

par(mfrow=c(3,3), oma=c(2,2,3,3), mar=c(4,4,1,0))
x=c(-4,-0.5)
y=c(0,1100)
bk=40

a<-read.table("D-High.DS_only", head=F)
c<-read.table("D-High.LFMD_only", head=F)
b<-read.table("D-High.overlap", head=F)

hist(log10(a[,1]), breaks=bk, xlab=NULL, ylab="Frequency (D-High)", main="DS+/LFMD-",xlim=x, ylim=y)
hist(log10(b[,1]), breaks=bk, xlab=NULL, ylab=NULL, main="DS+/LFMD+",xlim=x, ylim=y)
hist(log10(c[,1]), breaks=bk, xlab=NULL, ylab=NULL, main="DS-/LFMD+",xlim=x, ylim=y)

a<-read.table("ND-High.DS_only", head=F)
c<-read.table("ND-High.LFMD_only", head=F)
b<-read.table("ND-High.overlap", head=F)


hist(log10(a[,1]), breaks=bk, xlab=NULL, ylab="Frequency (ND-High)", main=NULL,xlim=x, ylim=y)
hist(log10(b[,1]), breaks=bk, xlab=NULL, ylab=NULL, main=NULL,xlim=x, ylim=y)
hist(log10(c[,1]), breaks=bk, xlab=NULL, ylab=NULL, main=NULL,xlim=x, ylim=y, col="red")

a<-read.table("ND-Low.DS_only", head=F)
c<-read.table("ND-Low.LFMD_only", head=F)
b<-read.table("ND-Low.overlap", head=F)

hist(log10(a[,1]), breaks=bk, xlab=NULL, ylab="Frequency (ND-Low)", main=NULL,xlim=x, ylim=y)
hist(log10(b[,1]), breaks=bk, xlab="log10(AAF)", ylab=NULL, main=NULL,xlim=x, ylim=y)
hist(log10(c[,1]), breaks=bk, xlab=NULL, ylab=NULL, main=NULL,xlim=x, ylim=y)

dev.off()
