pdf(file="CG2AT.ND-High.pdf", width=10, height=6)

par(mfrow=c(1,3), oma=c(2,2,3,3), mar=c(4,4,1,0))
x=c(-4,-0.5)
y=c(0,50)
bk=20

a<-read.table("ND-High.DS_only.mut.CG2AT", head=F)
c<-read.table("ND-High.LFMD_only.mut.CG2AT", head=F)
b<-read.table("ND-High.overlap.mut.CG2AT", head=F)

hist(log10(a[,6]), breaks=15, xlab=NULL, ylab="Frequency (D-High, C:G>A:T)", main="DS+/LFMD-",xlim=x, ylim=y)
hist(log10(b[,6]), breaks=30, xlab="log10(AAF)", ylab=NULL, main="DS+/LFMD+",xlim=x, ylim=y)
hist(log10(c[,6]), breaks=20, xlab=NULL, ylab=NULL, main="DS-/LFMD+",xlim=x, ylim=y, col="red")

dev.off()
