library(ggplot2)
library(svglite)
library(dplyr)
library(reticulate)
library(lattice)
library(gridExtra)
library(colormap)
######## Heatmaps ########
# cell count tables are pre-reordered for visualization purposes
# load tables and convert to proportions
cells_csigep = read.table("benchmarking/summary/cell_counts_CSI_GEP.txt", sep = "\t", header = 1, row.names = 1)
cells_csigep = as.matrix(cells_csigep)
colnames(cells_csigep) = c(1:40)
cells_csigep_proportion = cells_csigep/max(cells_csigep)

cells_scvi = read.table("benchmarking/summary/cell_counts_scvi.txt", sep = "\t", header = 1, row.names = 1)
cells_scvi = as.matrix(cells_scvi)
colnames(cells_scvi) = c(1:40)
cells_scvi_proportion = cells_scvi/max(cells_scvi)

cells_sagpu = read.table("benchmarking/summary/cell_counts_sagpu.txt", sep = "\t", header = 1, row.names = 1)
cells_sagpu = as.matrix(cells_sagpu)
colnames(cells_sagpu) = c(1:40)
cells_sagpu_proportion = cells_sagpu/max(cells_sagpu)

cells_inmf = read.table("benchmarking/summary/cell_counts_inmf.txt", sep = "\t", header = 1, row.names = 1)
cells_inmf = as.matrix(cells_inmf)
colnames(cells_inmf) = c(1:40)
cells_inmf_proportion = cells_inmf/max(cells_inmf)

numberOfBreaks <- 6
brksUniv <- seq(0, 1, length.out=numberOfBreaks)
myColorkey <- list(at=brksUniv, labels=list(at=brksUniv, labels=round(brksUniv,1)))
p1 <- levelplot(cells_csigep_proportion, margin=FALSE, xlab=NULL, ylab=NULL, aspect = "fill",
	main=list('CSI-GEP',side=0.5, cex=0.8, line=0.5), scales=list(draw=FALSE), at=brksUniv, colorkey=FALSE)
p2 <- levelplot(cells_inmf_proportion, margin=FALSE, xlab=NULL, ylab=NULL, aspect = "fill",
	main=list('iNMF',side=0.5, cex=0.8, line=0.5), scales=list(draw=FALSE), at=brksUniv, colorkey=FALSE)
p3 <- levelplot(cells_sagpu_proportion, margin=FALSE, xlab=NULL, ylab=NULL, aspect = "fill",
	main=list('SA-GPU',side=0.5, cex=0.8, line=0.5), scales=list(draw=FALSE), at=brksUniv, colorkey=FALSE)
p4 <- levelplot(cells_scvi_proportion, margin=FALSE, xlab=NULL, ylab=NULL, aspect = "fill",
	main=list('scVI',side=0.5, cex=0.8, line=0.5), scales=list(draw=FALSE), at=brksUniv, colorkey=FALSE)


color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
	scale = (length(lut)-1)/(max-min)
	par(mar=c(2,1,0,1))
	plot(c(min,max), c(-1,1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
	axis(1, ticks, las=1, tick=FALSE)
	for (i in 1:(length(lut)-1)) {
		y = (i-1)/scale + min
		rect(y,-2,y+1/scale, -1., col=lut[i], border=NA)
	}	
}
colors = colormap(colormaps[[12]]) # "YIGnBu"

pdf("benchmarking/summary/heatmaps.pdf", width = 16, height = 8) 
grid.arrange(p1, p2, p4, p3, color.bar(colors,1,0,nticks=6), ncol=2, nrow = 3, 
              layout_matrix = rbind(c(1,2), c(3,4), c(5,5)),
              widths = c(2.1, 2.6), heights = c(1.5, 1.5, 0.2))
dev.off()

######## Customized metrics ##########
N_scvi = dim(cells_csigep)[1]
N_csigep = dim(cells_csigep)[1]
N_sagpu = dim(cells_sagpu)[1]
N_inmf = dim(cells_inmf)[1]

# noise if the max % < 0.5
Nnoise_scvi = sum(apply(cells_scvi_proportion, 1, max)<0.5)
Nnoise_csigep = sum(apply(cells_csigep_proportion, 1, max)<0.5)
Nnoise_sagpu = sum(apply(cells_sagpu_proportion, 1, max)<0.5)
Nnoise_inmf = sum(apply(cells_inmf_proportion, 1, max)<0.5)

# falsely split if matching more than one latent representations
Nsplit_scvi = sum(apply(cells_scvi_proportion,2, FUN=function(x) sum(x>0.5)) > 1)
Nsplit_csigep = sum(apply(cells_csigep_proportion,2, FUN=function(x) sum(x>0.5)) > 1)
Nsplit_sagpu = sum(apply(cells_sagpu_proportion,2, FUN=function(x) sum(x>0.5)) > 1)
Nsplit_inmf = sum(apply(cells_inmf_proportion,2, FUN=function(x) sum(x>0.5)) > 1)

# falsely combined if more than one cell types matching to a latent representations
Ncombine_scvi = sum(apply(cells_scvi_proportion,1, FUN=function(x) sum(x>0.5)) > 1)
Ncombine_csigep = sum(apply(cells_csigep_proportion,1, FUN=function(x) sum(x>0.5)) > 1)
Ncombine_sagpu = sum(apply(cells_sagpu_proportion,1, FUN=function(x) sum(x>0.5)) > 1) 
Ncombine_inmf = sum(apply(cells_inmf_proportion,1, FUN=function(x) sum(x>0.5)) > 1)

df = data.frame(Label = rep(c("Normal", "Noise", "False split", "False combine"),4), 
	Models = c(rep("CSI-GEP",4), rep("iNMF",4), rep("SA-GPU",4), rep("ScVI",4)),
	Counts = c(N_csigep-Nnoise_csigep-Nsplit_csigep-Ncombine_csigep,Nnoise_csigep,Nsplit_csigep,Ncombine_csigep,
		N_inmf-Nnoise_inmf-Nsplit_inmf-Ncombine_inmf,Nnoise_inmf,Nsplit_inmf,Ncombine_inmf,
		N_sagpu-Nnoise_sagpu-Nsplit_sagpu-Ncombine_sagpu,Nnoise_sagpu,Nsplit_sagpu,Ncombine_sagpu,
		N_scvi-Nnoise_scvi-Nsplit_scvi-Ncombine_scvi,Nnoise_scvi,Nsplit_scvi,Ncombine_scvi))

p5 = ggplot(df, aes(fill = Label, y=Counts, x=Models)) + 
    geom_bar(stat="identity")  +
    theme(axis.title.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    	panel.background = element_blank(), axis.line = element_line(colour = "black"),
    	legend.position = c(0.4,0.3), text = element_text(size = 16), legend.background=element_rect(fill = alpha("white", 0.5))) +
    labs(y = "Number of GEPs", x = "", fill = "Customized metrics")


######## Simplified silouette ########
# function: simplified silhouette with adjusted cosine similarity
pairwiseMean = function(centroid, expression){
  return(mean(cbind(centroid,expression),na.rm = TRUE))
}

SS_cell = function(cellname, centroids, ranks){
  m = apply(centroids, 1, pairwiseMean, counts[cellname,])
  similarity = rowSums(t(t(centroids-m)*as.vector(counts[cellname,]-m)))/norm(as.matrix(counts[cellname,]-m),type='2')/sqrt(rowSums((centroids-m)^2))
  dist = 1-similarity
  a = dist[ranks[cellname,"gep"]]
  b = min(dist[-ranks[cellname,"gep"]],na.rm=TRUE)
  SS = (b-a)/max(a,b)
 return(SS)
}


# load count matrix and sample 5,000 cells to calculate on 
np <- import("numpy")
counts <- np$load("simulated_data/counts.npz", allow_pickle = TRUE)
counts = counts[['data']]
rownames(counts) = paste0("Cell", 1:nrow(counts))
metaCells <- np$load("simulated_data/cellparams.npz", allow_pickle = TRUE)
metaCells = metaCells[['data']]
rownames(metaCells) = paste0("Cell", 1:nrow(metaCells))

ss_names = paste0("Cell", sample(1:nrow(counts),5000,replace=FALSE))


####### csigep ########
# load decomposed matrices
gepindex <- readRDS("CSI-GEP/output_for_benchmarking/GEP_results/Splatter_250k_gepindex.RDS")
usages_split0 <- read.csv("CSI-GEP/output_for_benchmarking/cNMF_split0/cNMF_split0.usages.k_47.dt_0_1.consensus.txt", 
                          header = TRUE, sep = "\t", row.names = 1)
usages_split0_combined <- matrix(nrow = dim(usages_split0)[1], ncol = 41)

for (i in 1:length(gepindex)){
  ls = gepindex[[i]]$Split0
  if (length(ls)>1){
    usages_split0_combined[,i] <- rowMeans(usages_split0[,ls])
  }
  else{
    usages_split0_combined[,i] <- usages_split0[,ls]
  }
}
rownames(usages_split0_combined) = rownames(usages_split0)

usages_split1 <- read.csv("CSI-GEP/output_for_benchmarking/cNMF_split1/cNMF_split1.usages.k_47.dt_0_1.consensus.txt", 
                          header = TRUE, sep = "\t", row.names = 1)
usages_split1_combined <- matrix(nrow = dim(usages_split1)[1], ncol = 41)

for (i in 1:length(gepindex)){
  ls = gepindex[[i]]$Split1
  if (length(ls)>1){
    usages_split1_combined[,i] <- rowMeans(usages_split1[,ls])
  }
  else{
    usages_split1_combined[,i] <- usages_split1[,ls]
  }
}
rownames(usages_split1_combined) = rownames(usages_split1)

usages_combined <- rbind(usages_split0_combined, usages_split1_combined)
usages_combined = usages_combined[row.names(metaCells), ] # reorder ranks
ranks_csigep <- data.frame(gep = apply(usages_combined, 1, which.max), score = apply(usages_combined, 1, max))

# compute centroids of each GEP
centroids_csigep = matrix(NA, nrow = ncol(usages_combined), ncol = 25000)
for (i in 1:ncol(usages_combined)){
  centroids_csigep[i,] = colMeans(counts[rownames(ranks_csigep[ranks_csigep$gep==i & ranks_csigep$score > mean(ranks_csigep$score[ranks_csigep$gep==i]),]),])
}

# call simplified silouette
SS_ss = sapply(ss_names, SS_cell, centroids_csigep, ranks_csigep)
ss_csigep = mean(SS_ss)



####### inmf #######
usages = read.table("benchmarking/inmf/output/iNMF_usages.k_52.txt", sep = " ")
rownames(usages) = rownames(metaCells) 
ranks_inmf <- data.frame(gep = apply(usages, 1, which.max), score = apply(usages, 1, max))
centroids_inmf = matrix(NA, nrow = ncol(usages), ncol = 25000)
for (i in 1:ncol(usages)){
  centroids_inmf[i,] = colMeans(counts[rownames(ranks_inmf[ranks_inmf$gep==i & ranks_inmf$score > mean(ranks_inmf$score[ranks_inmf$gep==i]),]),])
}

SS_ss = sapply(ss_names, SS_cell, centroids_inmf, ranks_inmf)
ss_inmf = mean(SS_ss)



####### sagpu #######
usages <- read.csv("benchmarking/sa_gpu/output/output_W.txt", header = TRUE, sep = "\t", row.names = 1)
rownames(usages) = rownames(metaCells) 
ranks_sagpu <- data.frame(gep = apply(usages, 1, which.max), score = apply(usages, 1, max))
centroids_sagpu = matrix(NA, nrow = ncol(usages), ncol = 25000)
for (i in 1:ncol(usages)){
  centroids_sagpu[i,] = colMeans(counts[rownames(ranks_sagpu[ranks_sagpu$gep==i & ranks_sagpu$score > mean(ranks_sagpu$score[ranks_sagpu$gep==i]),]),])
}

SS_ss = sapply(ss_names, SS_cell, centroids_sagpu, ranks_sagpu)
ss_sagpu = mean(SS_ss)


###### scvi ######
usages = read.table("benchmarking/scvi/output/scvi_representation_42.txt", sep = " ", header = 0)
rownames(usages) = rownames(metaCells) 
ranks_scvi <- data.frame(gep = apply(usages, 1, which.max), score = apply(usages, 1, max))
centroids_scvi = matrix(NA, nrow = ncol(usages), ncol = 25000)
for (i in 1:ncol(usages)){
  centroids_scvi[i,] = colMeans(counts[rownames(ranks_scvi[ranks_scvi$gep==i & ranks_scvi$score > mean(ranks_scvi$score[ranks_scvi$gep==i]),]),])
}

SS_ss = sapply(ss_names, SS_cell, centroids_scvi, ranks_scvi)
ss_scvi = mean(SS_ss)

ss_df = data.frame(Silhouette = c(ss_scvi, ss_sagpu, ss_inmf, ss_csigep), 
	model = c("ScVI", "SA-GPU", "iNMF","CSI-GEP" ))
ss_df$model = factor(ss_df$model, levels=c("CSI-GEP",  "iNMF", "SA-GPU","ScVI" ))

p6 = ggplot(ss_df, aes(x=as.factor(model), y = Silhouette, fill = as.factor(model) )) + 
  geom_bar(stat = 'identity') + geom_text(aes(label=round(Silhouette,2)), vjust= 1) +
  scale_fill_hue(c = 40) +
  theme(legend.position="none",
    text = element_text(size = 12),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y = "Simplified Silhouette Score")

######## Mutual information ###########
# function: mutual information based on contingency table
MI_contingency <- function(tbl)
{
  C1 = rowSums(tbl)
  C2 = colSums(tbl)
  N = sum(tbl)

   
  MI = matrix(NA, ncol = ncol(tbl), nrow = nrow(tbl))
  for (i in 1:length(C1)){
    for (j in 1:length(C2)){
      MI[i,j] = tbl[i,j]/N *(log2((tbl[i,j]+1)/(C1[i]*C2[j]+1)) + log2(N))
    }
  }
  
  return(sum(MI))
}


MI_csigep = MI_contingency(cells_csigep)
MI_inmf = MI_contingency(cells_inmf)
MI_sagpu = MI_contingency(cells_sagpu)
MI_scvi = MI_contingency(cells_scvi)

mi_df = data.frame(MI = c(MI_scvi, MI_sagpu, MI_inmf, MI_csigep), 
  model = c("ScVI", "SA-GPU", "iNMF","CSI-GEP" ))
mi_df$model = factor(mi_df$model, levels=c("CSI-GEP",  "iNMF", "SA-GPU","ScVI" ))

p7 = ggplot(mi_df, aes(x=as.factor(model), y = MI, fill = as.factor(model) )) + 
  geom_bar(stat = 'identity') + geom_text(aes(label=round(MI,2)), vjust=1) +
  scale_fill_hue(c = 40) +
  theme(legend.position="none",
    text = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y = "Adjusted Rand Index")




######## ARI ########
# function: adjusted Rand index based on contingency table
ARI_contingency <- function(tbl)
{
  C1 = rowSums(tbl)
  C2 = colSums(tbl)
  n = sum(tbl)
  Exp_RI = (sum(C1*(C1-1))*sum(C2*(C2-1))) / (2*n*(n-1))
  RI = sum(tbl*(tbl-1)/2)
  max_RI = (sum(C1*(C1-1))+sum(C2*(C2-1)))/4
  
  Adj_RI = (RI-Exp_RI)/(max_RI-Exp_RI)
  return(Adj_RI)
}

ARI_csigep = ARI_contingency(cells_csigep)
ARI_inmf = ARI_contingency(cells_inmf)
ARI_sagpu = ARI_contingency(cells_sagpu)
ARI_scvi = ARI_contingency(cells_scvi)

ari_df = data.frame(ARI = c(ARI_scvi, ARI_sagpu, ARI_inmf, ARI_csigep), 
  model = c("ScVI", "SA-GPU", "iNMF","CSI-GEP" ))
ari_df$model = factor(ari_df$model, levels=c("CSI-GEP",  "iNMF", "SA-GPU","ScVI" ))

p8 = ggplot(ari_df, aes(x=as.factor(model), y = ARI, fill = as.factor(model) )) + 
  geom_bar(stat = 'identity') + geom_text(aes(label=round(ARI,2)), vjust=1) + 
  scale_fill_hue(c = 40) +
  theme(legend.position="none",
    text = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y = "Adjusted Rand Index")


###########################################
pdf("benchmarking/summary/metrics.pdf") 
grid.arrange(p5,p6,p7,p8, ncol=2, nrow = 2, 
              layout_matrix = rbind(c(1,2),c(3,4)),
              widths = c(1,1), heights = c(1,1))
dev.off()
