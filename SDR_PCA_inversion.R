library(ggplot2)
library(ggpubr)
library(dplyr)

################# 887 wild-caught individuals by SDR #################
for (chr in c(12, 3)){
  title = ifelse(chr==3, "LG3:17260000-17340000", "LG12:1-16900000")
  data = paste0("Feng887_LG", chr, "sex_a1m3_LG", chr, "pop")
  pca = read.table(paste0(data, ".eigenvec"))
  miss = read.table(paste0(data, ".imiss"), header=1)
  names(pca) =c("ID", paste0("PC", 1:(ncol(pca)-1)))
  pca = merge(pca, miss, by.x="ID", by.y="INDV")
  pca = merge(pca, sif887, by.x="ID", by.y="SampleID")
  eig=read.table(paste0(data, ".eigenval"))
  
  if (chr==12) {
    assign(paste0("LG", chr, "sex_pca.plot"), ggplot() +
             geom_point(data=pca[pca$Population != "POL-GDY",], 
                        aes(x=PC1, y=PC2, fill=sex), size=2, stroke=0.3, alpha=0.5, shape=21) +
             geom_text(data=pca[pca$Population == "POL-GDY" & pca$sex_region != "LG12",], 
                       aes(x=PC1 + 0.005, y=PC2, color=sex, label=Population), size=2) +
             geom_text(data=pca[pca$Population == "POL-GDY" & pca$sex_region == "LG12",], 
                       aes(x=PC1, y=PC2, color=sex, label=ID), size=2) +
             scale_fill_manual(values = c("Female"="indianred", "Male"="royalblue")) +
             scale_color_manual(values = c("Female"="indianred", "Male"="royalblue"), guide="none") +
             geom_hline(yintercept=0, color="grey60", linewidth=0.1) + 
             geom_vline(xintercept=0, color="grey60", linewidth=0.1) +
             labs(x=paste0("PC1: ", round(eig[1,1]/sum(eig)*100), "%"), 
                  y=paste0("PC2: ", round(eig[2,1]/sum(eig)*100),"%"),
                  title=title) + theme_bw() +
             theme(legend.background = element_rect(fill = "transparent")))
  }
  if (chr==3) {
    assign(paste0("LG", chr, "sex_pca.plot"), ggplot() +
             geom_point(data=pca[pca$Population != "POL-GDY",], 
                        aes(x=PC1, y=PC2, fill=sex), size=2, stroke=0.3, alpha=0.5, shape=21) +
             geom_text(data=pca[pca$Population == "POL-GDY" & pca$sex == "Male" &
                                  pca$ID != "POL-GDY-16",], 
                       aes(x=PC1, y=PC2, color=sex, label=Population), size=2) +
             geom_text(data=pca[pca$Population == "POL-GDY" & pca$sex == "Female",], 
                       aes(x=PC1 - 0.02, y=PC2, color=sex, label=ID), size=2) +
             geom_text(data=pca[pca$ID == "POL-GDY-16",], 
                       aes(x=PC1 - 0.02, y=PC2, color=sex, label=ID), size=2) +
             scale_fill_manual(values = c("Female"="indianred", "Male"="royalblue")) +
             scale_color_manual(values = c("Female"="indianred", "Male"="royalblue"), guide="none") +
             geom_hline(yintercept=0, color="grey60", linewidth=0.1) + 
             geom_vline(xintercept=0, color="grey60", linewidth=0.1) +
             labs(x=paste0("PC1: ", round(eig[1,1]/sum(eig)*100), "%"), 
                  y=paste0("PC2: ", round(eig[2,1]/sum(eig)*100),"%"),
                  title=title) + theme_bw() +
             theme(legend.background = element_rect(fill = "transparent")))
  }
}

chr=3  
my3pop_colo = c("FRA-VEY*" = "#fdae61",
                "GBR-GRO*" = "#54aa47",
                "SCO-HAR*" = "#2b83ba")
tree.ids = c("FRA-VEY-1", "FRA-VEY-10", 
             "SCO-HAR-1", "SCO-HAR-11", "SCO-HAR-2", 
             "GBR-GRO-10", "GBR-GRO-11", "GBR-GRO-4")
my3pop = ggplot() +
  geom_point(data=pca[!(pca$ID %in% tree.ids), ], 
             aes(x=PC1, y=PC2), size=1.5, stroke=1, alpha=0.5, color="grey50") +
  geom_point(data=pca[pca$Population %in% c("FRA-VEY*", "GBR-GRO*", "SCO-HAR*") & !(pca$ID %in% tree.ids), ], 
             aes(x=PC1, y=PC2, color=Population), size=1.5, stroke=1, alpha=0.5) +
  geom_point(data=pca[pca$ID %in% tree.ids, ], 
             aes(x=PC1, y=PC2, fill=Population), shape=21, color="black", size=1.5, stroke=1, alpha=1) +
  scale_color_manual(values = my3pop_colo) +
  geom_hline(yintercept=0, color="grey60", linewidth=0.1) + 
  geom_vline(xintercept=0, color="grey60", linewidth=0.1) +
  labs(x=paste0("PC1: ", round(eig[1,1]/sum(eig)*100), "%"), 
       y=paste0("PC2: ", round(eig[2,1]/sum(eig)*100),"%"),
       title=title) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "transparent"))

### additional evidence for an inversion
het = read.table("Feng887_LG3sex_a1m3_LG3pop.het", header=1)
lg3_data = merge(pca, het, by.x = "ID", by.y = "INDV")
lg3_data$O_HET = (lg3_data$N_SITES - lg3_data$O.HOM.) / lg3_data$N_SITES

mean(lg3_data[lg3_data$sex == "Male" & lg3_data$Population != "POL-GDY", "PC1"])
#[1] -0.0589662
## make heterozygotes around zero
lg3_data$PC1 = lg3_data$PC1 + 0.0589662
## fold
lg3_data$PC1 = abs(lg3_data$PC1)

## remove spurious individuals 
p2.data = lg3_data[!(lg3_data$Population %in% c("SCO-HAR*", "GBR-GRO*", "FRA-VEY*")) & 
                     lg3_data$ID != "POL-GDY-16", ]
cor.test(x=p2.data$PC1, y=p2.data$O_HET, method = "pearson", alternative = "two.sided")
# 	Pearson's product-moment correlation
#data:  p2.data$PC1 and p2.data$O_HET
#t = -54.722, df = 163, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.9807120 -0.9645804
#sample estimates:
#  cor 
#-0.9738463 

p2.2 = ggplot(p2.data[order(p2.data$sex),], 
              aes(x=PC1, y=O_HET)) + 
  geom_point(aes(color=sex), alpha=0.6) +
  scale_color_manual(values=c("Female"="indianred", "Male"="royalblue")) +
  geom_smooth(method='lm', color="grey30", lwd=0.5, alpha=0.4) +
  geom_text(x=0.1, y=0.55, label="r=-0.97\np-value < 2.2e-16") +
  labs(y="observed heterozygosity", title = "",
       x="PC1 (centered and folded)")

p3 = ggplot() + 
  geom_point(data = lg3_data, aes(x=PC1, y=O_HET), color="grey60", alpha=0.8) +
  geom_point(data = lg3_data[lg3_data$Population %in% c("SCO-HAR*", "GBR-GRO*", "FRA-VEY*"),],aes(x=PC1, y=O_HET, color=Population), alpha=0.6) +
  scale_color_manual(values = my3pop_colo) +
  labs(y="observed heterozygosity", title="",
       x="PC1 (centered and folded)")

### Figure panels 
row1 = ggarrange(ncol=3, nrow=1, 
                 LG12sex_pca.plot + geom_text(aes(x=-0.02, y=0), label="XX", size=4) +
                   geom_text(aes(x=0.03, y=0), label="XY", size=4) + 
                   theme(legend.position = "none", title = element_text(size=10)), 
                 LG3sex_pca.plot + geom_text(aes(x=0.01, y=0.1), label="XX", size=4) +
                   geom_text(aes(x=-0.05, y=0.1), label="XY", size=4) + 
                   geom_text(aes(x=-0.17, y=0.04), label="Y'Y'?", size=4) + 
                   theme(legend.position = "none", title = element_text(size=10)),
                 my3pop + theme(legend.position= "none", title = element_text(size=10)),
                 labels = c("A)", "B)", "C)"))
row2 = ggarrange(ncol=2, nrow=1, 
                 p3 + theme(title = element_text(size=10)),
                 p2.2 + theme(title = element_text(size=10)), 
                 widths = c(1, 1), labels = c("D)", "E)"))
pdf("../../LG3_introgression_Figures/Fig_Feng887_pcasex_Y.pdf", width = 11, height = 6)
png("../../LG3_introgression_Figures/Fig_Feng887_pcasex_Y.png", width = 8, height = 5, units = "in", res=600)
ggarrange(ncol=1, nrow=2, row1, row2)
dev.off()
