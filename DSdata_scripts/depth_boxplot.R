pdf(file="depth_boxplot.pdf", width=16, height=10)


a<-read.table(file="depth_form.txt", head=T)
library(ggplot2)
library(ggpubr)

the <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
theme <- theme( text = element_text(size=16) )

compare <- list( c("ND-Low", "ND-High"), c("ND-Low","D-High"), c("ND-High","D-High") )

sig <- stat_compare_means(comparisons = compare, method="t.test", label="p.signif")
labs <- labs(y="Depth", x="Position on ChrM" )
mlabs <- labs(y="Depth", x="Position on ChrM" )

bo <- geom_boxplot(aes(group = cut_width(Position, 300)) )
dot <- geom_dotplot(aes(group = cut_width(Position, 300)) )

d <- ggplot(a[a$Label=="D-High", ], aes(Position, Depth, fill=Label))
m <- ggplot(a[a$Label=="ND-High", ], aes(Position, Depth, fill=Label))
n <- ggplot(a[a$Label=="ND-Low", ], aes(Position, Depth, fill=Label))

d + bo + ylim(0,250000) +theme
m + bo + ylim(0,250000) +theme
n + bo + ylim(0,250000) +theme

ggplot(a, aes(Class, Depth, fill=Class) ) + geom_violin() + theme

dev.off()

