pdf(file="data_stat.pdf", width=16, height=10)


a<-read.table(file="data_stat_form.txt", head=T)
library(ggplot2)
library(ggpubr)

the <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
theme <- theme( text = element_text(size=16) )

compare <- list( c("ND-Low", "ND-High"), c("ND-Low","D-High"), c("ND-High","D-High") )
sig <- stat_compare_means(comparisons = compare, method="t.test", label="p.signif")

ggplot(a, aes(Label, Nread, fill=Class) ) + geom_violin() + theme 
ggplot(a, aes(Label, Nsscs, fill=Class) ) + geom_violin() + theme 
ggplot(a, aes(Label, Ndcs, fill=Class) ) + geom_violin() + theme
ggplot(a, aes(Label, Ratio, fill=Class) ) + geom_violin() + theme 

dev.off()

