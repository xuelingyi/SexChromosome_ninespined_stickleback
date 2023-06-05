library(ggplot2)
library(ggpubr)
library(scatterpie)
library(ggrepel)

ind237 = read.csv("marine237.csv")
pop12 = read.csv("marine_pop12.csv")
col.W = "#2a1657"
col.E = "#e9a64c"

##################################################################################################
############################################# pie charts on map ##########################################
source("../my.map.R")

## pie charts of SDR in males
pop12$male_LG3 = sapply(pop12$Population, function(x){ return(nrow(subset(ind237, Population == x & sex == "M" & sex_region == "LG3"))) })
pop12$male_LG12 = sapply(pop12$Population, function(x){ return(nrow(subset(ind237, Population == x & sex == "M" & sex_region == "LG12"))) })
my.map(pop12) + geom_scatterpie(data=pop12, aes(x=long_ling, y=lat_ling), cols = c("male_LG12", "male_LG3"), color=NA, alpha=0.8) + scale_fill_manual(values=c("male_LG3"=col.W, "male_LG12"=col.E))

## pie charts of sexes
pop12$male = sapply(pop12$Population, function(x){ return(nrow(subset(ind237, Population == x & sex == "M"))) })
pop12$female = sapply(pop12$Population, function(x){ return(nrow(subset(ind237, Population == x & sex == "F"))) })
my.map(pop12) + geom_scatterpie(data=pop12, aes(x=long_ling, y=lat_ling), cols = c("female", "male"), color=NA, alpha=0.8) + scale_fill_manual(values=c("male"="royalblue", "female"="indianred"))
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
############################################ plot admixture CV errors ################################################
datasets = c("marine237_autosome", "marine_male119_autosome", "marine_female118_autosome", 
             "marine_female118_LG12SDR", "marine_female118_LG3SDR")
files = c("./auto/adm/cv_k1.13.txt", "./auto/adm_m/adm_m_cv.txt", "./auto/adm_f/adm_f_cv.txt", 
          "./sex_region/adm/LG12sexfemale/LG12sexfemale_CV.txt", "./sex_region/adm/LG3sexfemale/LG3sexfemale_CV.txt")
for (i in 1:5){
  error = read.table(files[i])
  if(i == 1){ error = error[, 5:6] } else 
    { error = error[, 3:4] }
  colnames(error) = c("K", "CV_error")
  error$K = gsub("(K=", "", error$K, fixed=T)
  error$K = gsub("):", "", error$K, fixed=T)
  assign(paste0(datasets[i], "_CV.plot"), ggplot(data=error, aes(x=K, y=CV_error)) + geom_point(size=0.8) + stat_summary(fun=mean, geom="line", aes(group=1)) + scale_x_discrete(limits=factor(1:max(as.numeric(error$K)))) + labs(title=datasets[i], y="Cross-Validation error"))
}
ggarrange(nrow=2, ncol=3, labels = "AUTO",
          marine237_autosome_CV.plot, marine_male119_autosome_CV.plot, marine_female118_autosome_CV.plot, marine_female118_LG12SDR_CV.plot, marine_female118_LG3SDR_CV.plot)
##################################################################################################
######################################### ADMIXTURE plots #################################################

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


