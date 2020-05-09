pdf(file="adj_mutSpectrum.pdf", width=12, height=10)

a<-read.table(file="adj_mutSpectrum.txt", head=T)
library(ggplot2)
library(ggpubr)

#The unpaired two-samples Wilcoxon test (also known as Wilcoxon rank sum test or Mann-Whitney test) is a non-parametric alternative to the unpaired two-samples t-test, which can be used to compare two independent groups of samples. It’s used when your data are not normally distributed.

theme <- theme( text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))

compare <- list( c("ND-Low", "D-High"), c("ND-High","ND-Low" ), c("ND-High","D-High") )
sig <- stat_compare_means(comparisons = compare, method="wilcox.test", label="p.signif", hide.ns = TRUE,)
#sig <- stat_compare_means(comparisons = compare, method="t.test", label="p.signif", hide.ns = TRUE,)
mlabs <- labs(y="Adjusted # of mutations" )

a$Label <- factor(a$Label, levels=c("ND-Low", "ND-High", "D-High"))
a$Class <- factor(a$Class, levels=c("ND-Low_Hippocampus", "ND-Low_Parietal", "ND-High_Hippocampus", "ND-High_Parietal", "D-High_Hippocampus", "D-High_Parietal"))

# Class
ggplot(a[a$Tools=="LFMD", ], aes(Mutation, Count, fill=Class)) + geom_boxplot() + mlabs + labs(title="LFMD") + theme

# Label
ggplot(a[a$Tools=="LFMD", ], aes(Mutation, Count, fill=Label)) + geom_boxplot() + mlabs + labs(title="LFMD") + theme

ggboxplot(a, x="Label", y="Count", facet.by = "Mutation", fill="Label" ) +  mlabs + sig + theme
#ggboxplot(a[a$Sample != "558H",], x="Label", y="Count", facet.by = "Mutation", fill="Label" ) +  mlabs + sig + theme

dev.off()

