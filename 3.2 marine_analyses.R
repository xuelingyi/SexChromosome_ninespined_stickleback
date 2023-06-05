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
########################################### HIest ###########################################
library(HIest)
library(vcfR)
## all populations
endpops = "DENRUS"
f=10
## three populations
endpops = "GERSWE"
f=4

## run below for each dataset above
P = read.csv(paste0(endpops, "_freq.csv"), sep=" ")
p = P[P$abs.freq >= f/10, 1:4]
names(p) = c("Locus", "Allele",  "P1", "P2")
data = read.vcfR(paste0("marine237_a2m7_no3.12_LD2_a5g1_", endpops, "g0freq", f, ".recode.vcf")) 
G = t(data@gt)
G = G[-1, ]
if(!(dim(G)[2] == nrow(p))){print("P and G data do not match in nsnps!")}
nsnps = dim(G)[2] 
ninds = dim(G)[1] 
g = as.data.frame(NULL)
for(i in 1:nrow(G)){
    count = sapply(G[i,], FUN = function(x){
      if(!is.na(x)) {
        gt = unlist(strsplit(x, split="/"))
        x = as.numeric(gt[1]) + as.numeric(gt[2])}
      else{x = "NA"}
    })
    g = rbind(g, count)
}
row.names(g) = row.names(G)
colnames(g) = (data@fix)[,3] 
  
HI.ac <- HIest(g, p,  type="allele.count",
               iterations=500, surf=TRUE, startgrid=50)
HI.ac$SampleID = row.names(G)
HI.ac = merge(HI.ac, 
              ind237[, c("SampleID", "Species", "Population", "sex_region", "sex",
                           "adm2", "adm3", "adm4", "adm5")], 
              by="SampleID")
save(f, nsnps, ninds, p, g, HI.ac, 
     file=paste0(endpops, "g0freq", f, ".RData"))
ggplot(data=HI.ac) + 
    geom_segment(aes(x = 0, y = 0, xend = 0.5, yend = 1)) + 
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0)) + 
    geom_segment(aes(x = 0.5, y = 1, xend = 1, yend = 0)) +
    geom_point(aes(x=S, y=H, color=Population)) + 
    labs(x="S", y="HI", 
         title = paste0(endpops, "_freq", f/10, ": n=", ninds, ", snps=", nsnps)) +
    theme_classic()
}


