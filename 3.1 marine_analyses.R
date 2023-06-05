ind237 = read.csv("marine237.csv")
pop12 = read.csv("marine_pop12.csv")
library(ggplot2)
library(ggpubr)
library(scatterpie)
library(ggrepel)

##################################################################################################
############################################ plot kinship ###############################################
king = read.table("marine237_a2m7_no3.12_LD2_a5g1.kin0", header = F)
names(king) = c("FID1",	"IID1",	"FID2",	"IID2",	"NSNP",	"HETHET", "IBS0", "KINSHIP")
for (pop in unique(ind237$Population)) {
  assign(paste0(gsub("-", "_", pop), "_king"), ggplot() + 
           geom_histogram(data=king[king$IID1 %in% ind237[ind237$Population == pop, "SampleID"] & king$IID2 %in% ind237[ind237$Population == pop, "SampleID"], ],
                          aes(x=KINSHIP), fill="grey20") +
           labs(title=paste0(pop, ", N=", nrow(ind237[ind237$Population == pop, ]))) +
           theme_bw() + theme(plot.title=element_text(size=8), axis.title = element_text(size=6), axis.text = element_text(size=6)) 
           )
}
ggarrange(nrow=3, ncol=4, 
          RUS_LEV_king, FIN_KIV_king, FIN_HEL_king, FIN_HAM_king, 
          FIN_SEI_king, FIN_TVA_king, SWE_GOT_king, SWE_FIS_king, 
          SWE_BOL_king, POL_GDY_king, GER_RUE_king, DEN_NOR_king)
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
