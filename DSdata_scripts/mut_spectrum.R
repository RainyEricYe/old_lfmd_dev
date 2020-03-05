pdf(file="mut_spectrum.pdf", width=12, height=8)

a<-read.table(file="mut_info.txt", head=T)
library(ggplot2)
library(ggpubr)


theme <- theme( text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))

compare <- list( c("DS", "LFMD"), c("UMI","LFMD"), c("UniC","LFMD") )
sig <- stat_compare_means(comparisons = compare, method="t.test", label="p.signif")
labs <- labs(y="# of mutations", x="Tools" )
mlabs <- labs(y="# of mutations", x="Mutation" )

# Class
ggplot(a[a$Tools=="LFMD", ], aes(Mutation, Count, fill=Class)) + geom_boxplot() + mlabs + labs(title="LFMD") + theme
ggplot(a[a$Tools=="DS", ], aes(Mutation, Count, fill=Class)) + geom_boxplot() + mlabs + labs(title="DS") + theme

# Label
ggplot(a[a$Tools=="LFMD", ], aes(Mutation, Count, fill=Label)) + geom_boxplot() + mlabs + labs(title="LFMD") + theme
ggplot(a[a$Tools=="DS", ], aes(Mutation, Count, fill=Label)) + geom_boxplot() + mlabs + labs(title="DS") + theme

# only CG>AT
ggplot(a[ (a$Dementia=="ND" & a$Mutation=="C:G>A:T"), ], aes(Tools, Count, fill=Class )) + geom_boxplot() + labs + labs(title="C:G>A:T") + theme 

ggboxplot(a[ a$Mutation=="C:G>A:T", ], x="Tools", y="Count", facet.by = "Label", fill="Tools" ) +  labs + sig +labs(title="C:G>A:T") + theme

ggboxplot(a[ a$Mutation=="C:G>A:T", ], x="Tools", y="Count", facet.by = "Class", fill="Tools" ) +  labs + sig + labs(title="C:G>A:T") + theme


dev.off()

