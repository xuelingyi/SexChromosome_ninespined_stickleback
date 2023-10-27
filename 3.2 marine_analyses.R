library(ggplot2)
library(ggpubr)
library(scatterpie)
library(ggrepel)

ind237 = read.csv("marine237.csv")
pop12 = read.csv("marine_pop12.csv")
col.W = "#2a1657"
col.E = "#e9a64c"
col.Sweden.west = "#665393"
col.Germany="#b8add1"
col.East2="#ffc475"
col.Russia="#a66610"
col.7 = "#fff7bc"
col.pop = c("RUS-LEV"="#6a3d9a", 
            "FIN-HAM"="#cab2d6", "FIN-KIV"="#a6cee3", "FIN-HEL"="#1f78b4", "FIN-TVA"="#33a02c", "FIN-SEI"="#8dd3c7",
            "SWE-BOL"="#b2df8a", "SWE-GOT"="#ffff99","POL-GDY"="#fdbf6f", "GER-RUE"="#ff7f00", "SWE-FIS"="#fb9a99", 
            "DEN-NOR"= "#e31a1c")
order.adm=read.table("order.adm.txt")
order.adm = order.adm$V1
##################################################################################################
############################################# pie charts on map ##########################################
## pie charts of SDR in males
pop12$male_LG3 = sapply(pop12$Population, function(x){ return(nrow(subset(ind237, Population == x & sex == "M" & sex_region == "LG3"))) })
pop12$male_LG12 = sapply(pop12$Population, function(x){ return(nrow(subset(ind237, Population == x & sex == "M" & sex_region == "LG12"))) })
ggplot() + geom_scatterpie(data=pop12, aes(x=long_ling, y=lat_ling), cols = c("male_LG12", "male_LG3"), color=NA, alpha=0.8) + scale_fill_manual(values=c("male_LG3"=col.W, "male_LG12"=col.E))

## pie charts of sexes
pop12$male = sapply(pop12$Population, function(x){ return(nrow(subset(ind237, Population == x & sex == "M"))) })
pop12$female = sapply(pop12$Population, function(x){ return(nrow(subset(ind237, Population == x & sex == "F"))) })
ggplot() + geom_scatterpie(data=pop12, aes(x=long_ling, y=lat_ling), cols = c("female", "male"), color=NA, alpha=0.8) + scale_fill_manual(values=c("male"="royalblue", "female"="indianred"))
##################################################################################################
############################################ plot ADMIXTURE results ################################################
datasets = c("marine237_autosome", "marine_male119_autosome", "marine_female118_autosome")

### CV errors
files = c("cv_k1.13.txt", "adm_m_cv.txt", "adm_f_cv.txt")
for (i in 1:3){
  error = read.table(files[i])
  if(i == 1){ error = error[, 5:6] } else 
    { error = error[, 3:4] }
  colnames(error) = c("K", "CV_error")
  error$K = gsub("(K=", "", error$K, fixed=T)
  error$K = gsub("):", "", error$K, fixed=T)
  assign(paste0(datasets[i], "_CV.plot"), ggplot(data=error, aes(x=K, y=CV_error)) + geom_point(size=0.8) + stat_summary(fun=mean, geom="line", aes(group=1)) + scale_x_discrete(limits=factor(1:max(as.numeric(error$K)))) + labs(title=datasets[i], y="Cross-Validation error"))
}
ggarrange(nrow=2, ncol=3, labels = "AUTO",
          marine237_autosome_CV.plot, marine_male119_autosome_CV.plot, marine_female118_autosome_CV.plot)

## admixture plots
library(reshape2)
paths = c("1675051151/K=", "1677129050/K=", "1677128605/K=") # CLUMPAK
IDlists = c("ID", "ID.m", "ID.f")

for (i in 1:3){
  ID = read.table(IDlists[i])
  max = ifelse(i==1, 7, 6)
  pdf(paste0(datasets[i], "_clum_plot.pdf"))
  for(K in 2:max){
    clu_adm = read.table(paste0(paths[i], K, "/CLUMPP.files/ClumppIndFile.output"), row.names = 1)
    colnames(clu_adm) = c("sample", "zero", "pop", "colon", paste0("cluster",1:K))
    clu_adm$sample = ID$V1
    clu_adm = merge(clu_adm, ind237[, c("SampleID", "Population")], by.x ="sample", by.y="SampleID")
    clu_adm_ordered = clu_adm
    clu_adm_ordered$sample = factor(clu_adm_ordered$sample, levels=order.adm)
    clu_adm_plot = clu_adm_ordered[, !colnames(clu_adm_ordered) %in% c("zero", "Population", "pop", "colon")]
    clu_adm_plot = melt(clu_adm_plot, id.vars="sample")
    p = ggplot(clu_adm_plot, aes(x=sample, y=value, fill=variable)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle=90, size=3), axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0)) + labs(title=paste0("K", K), y="Probability")
    assign(paste0(datasets[i], "_K", K, "_clu"), p)
    print(p)}
  dev.off()}

top.themes = theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

ggarrange(nrow=6, ncol=1, heights = c(1,1,1,1,1,1.5), 
          marine237_autosome_K2_clu + top.themes + scale_fill_manual(values=c(col.W, col.E)),
          marine237_autosome_K3_clu + top.themes + scale_fill_manual(values=c(col.Sweden.west, col.W, col.E)),
          marine237_autosome_K4_clu + top.themes + scale_fill_manual(values=c(col.W, col.E, col.Russia, col.Sweden.west)),
          marine237_autosome_K5_clu + top.themes + scale_fill_manual(values=c(col.Sweden.west, col.E, col.Russia, col.W, col.Germany)),
          marine237_autosome_K6_clu + top.themes + scale_fill_manual(values=c(col.Sweden.west, col.East2, col.Germany, col.E, col.Russia, col.W)),
          marine237_autosome_K7_clu + theme(legend.position = "none") + scale_fill_manual(values=c(col.Russia, col.W, col.E, col.7, col.Sweden.west, col.Germany, col.East2)))

a = ggarrange(nrow=5, ncol=1, heights = c(1,1,1,1,1.5), 
              marine_male119_autosome_K2_clu + top.themes + scale_fill_manual(values=c(col.E, col.W)),
              marine_male119_autosome_K3_clu + top.themes + scale_fill_manual(values=c(col.E, col.W, col.Russia)),
              marine_male119_autosome_K4_clu + top.themes + scale_fill_manual(values=c(col.Sweden.west, col.W, col.E, col.Russia)),
              marine_male119_autosome_K5_clu + top.themes + scale_fill_manual(values=c(col.Germany, col.W, col.E, col.Russia, col.Sweden.west)),
              marine_male119_autosome_K6_clu + theme(legend.position = "none") + scale_fill_manual(values=c(col.Sweden.west, col.Russia, col.Germany, col.W, col.East2, col.E)))

b = ggarrange(nrow=5, ncol=1, heights = c(1,1,1,1,1.5), 
              marine_female118_autosome_K2_clu + top.themes + scale_fill_manual(values=c(col.E, col.W)),
              marine_female118_autosome_K3_clu + top.themes + scale_fill_manual(values=c(col.E, col.W, col.Sweden.west)),
              marine_female118_autosome_K4_clu + top.themes + scale_fill_manual(values=c(col.E, col.Russia, col.Sweden.west, col.W)),
              marine_female118_autosome_K5_clu + top.themes + scale_fill_manual(values=c(col.Sweden.west, col.W, col.Germany, col.Russia, col.E)),
              marine_female118_autosome_K6_clu + theme(legend.position = "none") + scale_fill_manual(values=c(col.East2, col.Russia, col.E, col.Sweden.west, col.W, col.Germany)))

ggarrange(ncol = 2, nrow=1, labels = datasets[2:3], a, b)
##################################################################################################
########################################### PCA #############################################
datasets = c("marine237_a2m7_no3.12_LD2_a5g1", "marine237_a2m7_no3.12_LD2_a5g1_m", "marine237_a2m7_no3.12_LD2_a5g1_f")
titles = c("marine237_a2m7_no3.12_LD2_a5g1", "marine_119males_autosome", "marine_118females_autosome")

for (i in 1:3){
  title=titles[i]
  pca = read.table(paste0(datasets[i], ".eigenvec"))
  eig=read.table(paste0(datasets[i], ".eigenval"))
  if(i == 1) { miss = read.table(paste0(datasets[i], ".vcf.imiss"), header=1)}
  if(i>=4){ miss = read.table(paste0(datasets[i], ".imiss"), header=1)}
  if(i<4){ names(pca) =c("FID", "ID", paste0("PC", 1:(ncol(pca)-2))) }
  if(i>=4){ names(pca) =c("ID", paste0("PC", 1:(ncol(pca)-1))) }
  pca = merge(pca, miss, by.x="ID", by.y="INDV")
  pca = merge(pca, ind237[, c("SampleID", "Population", "sex_region", "sex", "adm5")], by.x="ID", by.y="SampleID")
  
  assign(paste0(datasets[i], "_pca_miss.plot"), ggplot(pca, aes(x=PC1, y=PC2, color=F_MISS, shape=sex)) +
           scale_shape_manual(values=c("F"=16, "M"=17))+
           geom_point(size=2.5, stroke=1, alpha=0.5) +
           geom_hline(yintercept=0, color="grey60", size=0.1) + 
           geom_vline(xintercept=0, color="grey60", size=0.1) +
           labs(x=paste0("PC1: ", round(eig[1,1]/sum(eig)*100), "%"), y=paste0("PC2: ", round(eig[2,1]/sum(eig)*100),"%"), title=title) + theme_bw() +
           theme(legend.background = element_rect(fill = "transparent")))
  if (i<4){
    assign(paste0(datasets[i], "_pca_pop.plot"), ggplot() +
             geom_point(data=pca, aes(x=PC1, y=PC2, color=Population, shape=adm5), size=2.5) +
             scale_color_manual(values = col.pop) +
             scale_shape_manual(values=c("East"=17, "Denmark"=16, "Germany"=15, "Russia"=3, "Sweden.west"=7))+
             geom_hline(yintercept=0, color="grey60", size=0.1) + 
             geom_vline(xintercept=0, color="grey60", size=0.1) +
             labs(x=paste0("PC1: ", round(eig[1,1]/sum(eig)*100), "%"), y=paste0("PC2: ", round(eig[2,1]/sum(eig)*100),"%"), title=title) + theme_bw() +
             theme(legend.background = element_rect(fill = "transparent"))) }
  if(i >=4){
    assign(paste0(datasets[i], "_pca_pop.plot"), ggplot() +
             geom_point(data=pca, aes(x=PC1, y=PC2, fill=Population, shape=sex_region), size=2.5) +
             scale_fill_manual(values = col.pop) +
             scale_shape_manual(values=c("LG3"=24, "LG12"=21))+
             geom_hline(yintercept=0, color="grey60", size=0.1) + 
             geom_vline(xintercept=0, color="grey60", size=0.1) +
             labs(x=paste0("PC1: ", round(eig[1,1]/sum(eig)*100), "%"), y=paste0("PC2: ", round(eig[2,1]/sum(eig)*100),"%"), title=title) + theme_bw() +
             theme(legend.background = element_rect(fill = "transparent"))) }
  }
ggarrange(marine237_a2m7_no3.12_LD2_a5g1_pca_miss.plot, marine237_a2m7_no3.12_LD2_a5g1_pca_pop.plot,
          marine237_a2m7_no3.12_LD2_a5g1_m_pca_pop.plot, marine237_a2m7_no3.12_LD2_a5g1_f_pca_pop.plot,
          labels = "AUTO")
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
    labs(x="S", y="HI", title = paste0(endpops, "_freq", f/10, ": n=", ninds, ", snps=", nsnps)) +
    theme_classic()
}
##################################################################################################
############################################### IBDseq metadata ###############################################
## geographic distance in km
library("geosphere") 
eucl = distm(pop12[, c("long_ling", "lat_ling")], pop12[, c("long_ling", "lat_ling")]) ##in meter
row.names(eucl) = pop12$Population
colnames(eucl) = pop12$Population

## create metadata
ibd_metadata = function(myibd, sif, eucl, ...){ ## read in ibd and add pop, sex, sex_region
  ibd = read.table(myibd)
  names(ibd) = c("ind1", "haplotype1", "ind2", "haplotype2", "chr", "start", "end", "LOD")
  ibd$IBD_length = ibd$end - ibd$start
  ibd$pop1 = NULL
  ibd$pop2 = NULL
  ibd$geography.km = NULL
  ibd$sex1 = NULL
  ibd$sex2 = NULL
  ibd$mt1 = NULL
  ibd$mt2 = NULL
  for (i in 1:nrow(ibd)) {
    ibd[i, "pop1"] = sif[sif$SampleID == ibd[i, "ind1"], "Population"]
    ibd[i, "pop2"] = sif[sif$SampleID == ibd[i, "ind2"], "Population"]
    ibd[i, "geography.km"] = eucl[ibd[i, "pop1"], ibd[i, "pop2"]]/1000 #change meter to km
    ibd[i, "sex1"] = sif[sif$SampleID == ibd[i, "ind1"], "sex"]
    ibd[i, "sex2"] = sif[sif$SampleID == ibd[i, "ind2"], "sex"]
  }
  return(ibd)
}
auto19.ibd = ibd_metadata("marine237_a2m7_no3.12_LD2_a5g1.lod4.auto.ibd", sif=ind237, eucl = eucl)
save(auto19.ibd, file="IBDseq.RData")

lg3.ibd = ibd_metadata("LG3.marine237.a2m7.lod4.ibd", sif=ind237, eucl=eucl)
save(lg3.ibd, file="LG3_IBDseq.RData")
lg12.ibd = ibd_metadata("LG12.marine237.a2m7.lod4.ibd", sif=ind237, eucl=eucl)
save(lg12.ibd, file="LG12_IBDseq.RData")

####################################### IBD by chr ##############################################
## type of males
male3 = ind237[ind237$sex == "M" & ind237$sex_region == "LG3", "SampleID"]
male12 = ind237[ind237$sex == "M" & ind237$sex_region == "LG12", "SampleID"]
col.ibd = c("Male3 x Male3"="#08306b", "Male3 x Female"="#6baed6", "Male12 x Male12"="#662506", "Male12 x Female"="#fe9929", "Male3 x Male12" = "#41ab5d", "Female x Female" = "#dd3497")

load("LG3_IBDseq.RData")
load("LG12_IBDseq.RData")
load("IBDseq.RData")

## LOD-log10length
auto = ggplot(auto19.ibd) + geom_point(aes(x=log10(IBD_length.bp), y=LOD), size=1) + labs(title=paste0("19 autosomes #track: ", nrow(auto19.ibd))) + theme_bw()
lod3 = ggplot(lg3.ibd) + geom_point(aes(x=log10(IBD_length), y=LOD), size=1) + labs(title=paste0("LG3 #track: ", nrow(lg3.ibd)), x="log10(IBD_length.bp)") + theme_bw()
lod12 = ggplot(lg12.ibd) + geom_point(aes(x=log10(IBD_length), y=LOD), size=1) + labs(title=paste0("LG12 #track: ", nrow(lg12.ibd)), x="log10(IBD_length.bp)") + theme_bw()
ggarrange(nrow=1, ncol=3, auto, lod3, lod12)

auto19.ibd$type = apply(auto19.ibd, 1, function(x){
  if(x[1] %in% male3 & x[3] %in% male3) {return("Male3 x Male3")}
  if(x[1] %in% male12 & x[3] %in% male12) {return("Male12 x Male12")}
  if(x[1] %in% male3 & x[3] %in% male12) {return("Male3 x Male12")}
  if(x[1] %in% male12 & x[3] %in% male3) {return("Male3 x Male12")}
  if(x[13] == "F" & x[14] == "F") {return("Female x Female")}
  if(x[1] %in% male12 & x[14] == "F") {return("Male12 x Female")}
  if(x[13] == "F" & x[3] %in% male12) {return("Male12 x Female")}
  if(x[1] %in% male3 & x[14] == "F") {return("Male3 x Female")}
  if(x[13] == "F" & x[3] %in% male3) {return("Male3 x Female")}
})
lg3.ibd$type = apply(lg3.ibd, 1, function(x){
    if(x[1] %in% male3 & x[3] %in% male3) {return("Male3 x Male3")}
    if(x[1] %in% male12 & x[3] %in% male12) {return("Male12 x Male12")}
    if(x[1] %in% male3 & x[3] %in% male12) {return("Male3 x Male12")}
    if(x[1] %in% male12 & x[3] %in% male3) {return("Male3 x Male12")}
    if(x[13] == "F" & x[14] == "F") {return("Female x Female")}
    if(x[1] %in% male12 & x[14] == "F") {return("Male12 x Female")}
    if(x[13] == "F" & x[3] %in% male12) {return("Male12 x Female")}
    if(x[1] %in% male3 & x[14] == "F") {return("Male3 x Female")}
    if(x[13] == "F" & x[3] %in% male3) {return("Male3 x Female")}
  })
lg12.ibd$type = apply(lg12.ibd, 1, function(x){
  if(x[1] %in% male3 & x[3] %in% male3) {return("Male3 x Male3")}
  if(x[1] %in% male12 & x[3] %in% male12) {return("Male12 x Male12")}
  if(x[1] %in% male3 & x[3] %in% male12) {return("Male3 x Male12")}
  if(x[1] %in% male12 & x[3] %in% male3) {return("Male3 x Male12")}
  if(x[13] == "F" & x[14] == "F") {return("Female x Female")}
  if(x[1] %in% male12 & x[14] == "F") {return("Male12 x Female")}
  if(x[13] == "F" & x[3] %in% male12) {return("Male12 x Female")}
  if(x[1] %in% male3 & x[14] == "F") {return("Male3 x Female")}
  if(x[13] == "F" & x[3] %in% male3) {return("Male3 x Female")}
})

## v7 centromeric regions
centromere = read.csv("Kivikoski2021_TableS1_centromere.csv")
for(chr in c(1:21)){
  if(chr %in% c(3, 12)){
    ibd = get(paste0("lg", chr, ".ibd"))
  } else {
    ibd = auto19.ibd[auto19.ibd$chr == paste0("LG", chr), ]
  }
  ibd = ibd[order(ibd$LOD),]
  ibd$row = 1:nrow(ibd)
  
  if(chr == 3) {
    sdr = c(17260000, 17340000)
    assign(paste0("LG", chr, ".ibd"), ibd)}
  if(chr == 12) {
    sdr = c(1, 16900000)
    assign(paste0("LG", chr, ".ibd"), ibd)}
  
  assign(paste0("ibd", chr, ".plot"), ggplot(data=ibd) + 
           geom_segment(aes(x = start, y = row, xend = end, yend = row, color=type)) +
           scale_color_manual(values=col.ibd, breaks=c("Male3 x Male3", "Male3 x Female",  "Male12 x Male12", "Male12 x Female", "Male3 x Male12", "Female x Female")) +
           geom_vline(xintercept = centromere[centromere$Linkage_Group == chr, "Start.bp"]) +
           geom_vline(xintercept = centromere[centromere$Linkage_Group == chr, "End.bp"]) +
           theme_classic() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
           labs(x=paste0("LG", chr, " position (bp)"), y="IBD-like track", title=paste0("#tracks per Mbp: ", round(1000000*nrow(ibd)/max(c(ibd$start, ibd$end)),2)))
)
  if(chr %in% c(3, 12)) {
    assign(paste0("ibd", chr, ".plot"), ggplot(data=ibd) + 
             geom_segment(aes(x = start, y = row, xend = end, yend = row, color=type)) +
             scale_color_manual(values=col.ibd, breaks=c("Male3 x Male3", "Male3 x Female", "Male12 x Male12", "Male12 x Female", "Male3 x Male12", "Female x Female")) +
             geom_vline(xintercept = centromere[centromere$Linkage_Group == chr, "Start.bp"]) +
             geom_vline(xintercept = centromere[centromere$Linkage_Group == chr, "End.bp"]) +
             geom_vline(xintercept = sdr[1], color="red", linetype="dashed") +
             geom_vline(xintercept = sdr[2], color="red", linetype="dashed") +
             theme_classic() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
             labs(x=paste0("LG", chr, " position (bp)"), y="IBD-like track", title=paste0("#tracks per Mbp: ", round(1000000*nrow(ibd)/max(c(ibd$start, ibd$end)),2)))
    ) 
  }
}
ggarrange(ncol=3, nrow=7, 
          ibd1.plot, ibd2.plot, ibd3.plot,
          ibd4.plot, ibd5.plot, ibd6.plot,
          ibd7.plot, ibd8.plot, ibd9.plot,
          ibd10.plot, ibd11.plot, ibd12.plot,
          ibd13.plot, ibd14.plot, ibd15.plot,
          ibd16.plot, ibd17.plot, ibd18.plot,
          ibd19.plot, ibd20.plot, ibd21.plot)

####################################################################################################
################################################# IBD heatmap ##########################################
order.adm5 = c("DEN-NOR", "SWE-FIS", "GER-RUE", "POL-GDY", "SWE-GOT", "FIN-HEL", "FIN-SEI", "FIN-TVA", "SWE-BOL", "FIN-HAM", "FIN-KIV", "RUS-LEV")
color.adm5 = c(col.W, col.Sweden.west, col.Germany, col.E, col.E, col.E, col.E, col.E, col.E, col.E, col.E, col.Russia)
load("IBDseq.RData")
names(auto19.ibd)
auto19.ibd$IBD_length.Mbp = auto19.ibd$IBD_length.bp / 1000000
auto19.ibd$pair = paste0(auto19.ibd$ind1, "_", auto19.ibd$ind2, "_",
                         auto19.ibd$sex1, "_", auto19.ibd$sex2, "_",
                         auto19.ibd$sex_region1, "_", auto19.ibd$sex_region2, "_",
                         auto19.ibd$pop1, "_", auto19.ibd$pop2)
auto19.heatmap = as.data.frame(unique(auto19.ibd$pair))
names(auto19.heatmap) = "pair"
auto19.heatmap$sum.IBD_length.Mbp = sapply(auto19.heatmap$pair, function(x){return(sum(auto19.ibd[auto19.ibd$pair == x, "IBD_length.Mbp"]))} )
auto19.heatmap$ind1 = sapply(auto19.heatmap$pair, function(x) unlist(strsplit(x, split="_"))[1])
auto19.heatmap$ind2 = sapply(auto19.heatmap$pair, function(x) unlist(strsplit(x, split="_"))[2])
auto19.heatmap$sex1 = sapply(auto19.heatmap$pair, function(x) unlist(strsplit(x, split="_"))[3])
auto19.heatmap$sex2 = sapply(auto19.heatmap$pair, function(x) unlist(strsplit(x, split="_"))[4])
auto19.heatmap$sex_region1 = sapply(auto19.heatmap$pair, function(x) unlist(strsplit(x, split="_"))[5])
auto19.heatmap$sex_region2 = sapply(auto19.heatmap$pair, function(x) unlist(strsplit(x, split="_"))[6])
auto19.heatmap$pop1 = sapply(auto19.heatmap$pair, function(x) unlist(strsplit(x, split="_"))[7])
auto19.heatmap$pop2 = sapply(auto19.heatmap$pair, function(x) unlist(strsplit(x, split="_"))[8])

auto19.heatmap2 = auto19.heatmap
auto19.heatmap2$ind1 = auto19.heatmap$ind2
auto19.heatmap2$sex1 = auto19.heatmap$sex2
auto19.heatmap2$sex_region1 = auto19.heatmap$sex_region2
auto19.heatmap2$pop1 = auto19.heatmap$pop2
auto19.heatmap2$ind2 = auto19.heatmap$ind1
auto19.heatmap2$sex2 = auto19.heatmap$sex1
auto19.heatmap2$sex_region2 = auto19.heatmap$sex_region1
auto19.heatmap2$pop2 = auto19.heatmap$pop1

auto19.heatmap.full = rbind(auto19.heatmap, auto19.heatmap2)
any(duplicated(auto19.heatmap.full))

## all IBD-tracks 
ggplot(data=auto19.heatmap.full, aes(ind1, ind2, fill= log10(sum.IBD_length.Mbp))) + theme_bw() + geom_tile() + 
  theme(axis.text.x = element_text(angle=90, size=0.5), axis.text.y = element_text(size=0.5), axis.ticks = element_blank()) +
  scale_x_discrete(limits=order.adm) + scale_y_discrete(limits=order.adm) +
  scale_fill_gradient(low = "#efedf5", high = "#3f007d") 

## separate by sex
ind237$SampleID = factor(ind237$SampleID, levels = order.adm)
ind237 = ind237[order(ind237$SampleID),]
                             
female = ggplot(data=subset(auto19.heatmap.full, sex1=="F" & sex2=="F"), aes(ind1, ind2, fill= log10(sum.IBD_length.Mbp))) + 
  theme_bw() + geom_tile() + theme(axis.text.x = element_text(angle=90, size=0.5), axis.text.y = element_text(size=0.5),axis.ticks = element_blank()) +
  scale_x_discrete(limits=ind237[ind237$sex == "F", "SampleID"]) + scale_y_discrete(limits=ind237[ind237$sex == "F", "SampleID"]) +
  scale_fill_gradient(low = "#fee0d2", high = "#67000d")
                             
male = ggplot(data=subset(auto19.heatmap.full, sex1=="M" & sex2=="M"), aes(ind1, ind2, fill= log10(sum.IBD_length.Mbp))) + 
  theme_bw() + geom_tile() + theme(axis.text.x = element_text(angle=90, size=0.5), axis.text.y = element_text(size=0.5), axis.ticks = element_blank()) + 
  scale_x_discrete(limits=ind237[ind237$sex == "M", "SampleID"]) + scale_y_discrete(limits=ind237[ind237$sex == "M", "SampleID"]) +
  scale_fill_gradient(low = "#d3e1ff", high = "#061539")

ggarrange(ncol=2, nrow = 1, male, female)

##########################################################################################
############################################## IBD by genetic cluster #############################################
### adm5
range(auto19.ibd$geography.km)
range(auto19.ibd$IBD_length.Mbp)
adm5 = c("Denmark", "Sweden.west", "Germany", "East", "Russia")
for (i in 1:5){
  ID1 = ind237[ind237$adm5 == adm5[i], "SampleID"]
  for (j in i:5){
    ID2 = ind237[ind237$adm5 == adm5[j], "SampleID"]
    data = rbind(subset(auto19.ibd, ind1 %in% ID1 & ind2 %in% ID2), subset(auto19.ibd, ind1 %in% ID2 & ind2 %in% ID1))
    if(i==j){data = subset(auto19.ibd, ind1 %in% ID1 & ind2 %in% ID1)}
    assign(paste0("plot", i, j), ggplot(data) + geom_point(aes(x=geography.km, y=IBD_length.Mbp)) + 
             scale_x_continuous(limits = c(0, 1850))+ scale_y_continuous(limits = c(0, 31)) + theme_bw())
} }
##################################################################################################
