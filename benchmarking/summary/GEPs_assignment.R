library(reticulate)
np <- import("numpy")
metaCells <- np$load("simulated_data/cellparams.npz", allow_pickle = TRUE)
metaCells = metaCells[['data']]
rownames(metaCells) = paste0("Cell", 1:nrow(metaCells))

########################################
############## CSI-GEP #################
########################################
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
# usages_combiend generated from visualize.R
usages_combined = usages_combined[row.names(metaCells), ] # reorder ranks
ranks <- data.frame(gep = apply(usages_combined, 1, which.max), score = apply(usages_combined, 1, max))
table(ranks$gep)


high_score_tbl = matrix(nrow=41,ncol=40)
for (k in 1:41){
  tbl = table(factor(unlist(metaCells[rownames(metaCells) %in% rownames(usages_combined[usages_combined[,k]>mean(ranks$score[ranks$gep==k]),]),1]), 
                     levels = 1:40)) 
  high_score_tbl[k,] = tbl
}
write.table(high_score_tbl, file = "benchmarking/summary/cell_counts_CSI_GEP.txt", sep = "\t")



########################################
############## iNMF #################
########################################
usages = read.table("benchmarking/inmf/output/iNMF_usages.k_52.txt", sep = " ")
rownames(usages) = rownames(metaCells) 
ranks <- data.frame(gep = apply(usages, 1, which.max), score = apply(usages, 1, max))

high_score_tbl = matrix(nrow=52,ncol=40)
for (k in 1:52){
  tbl = table(factor(unlist(metaCells[rownames(metaCells) %in% rownames(usages[usages[,k]>mean(ranks[ranks$gep==k,]$score),]),1]), 
                     levels = 1:40)) 
  high_score_tbl[k,] = tbl
}
write.table(high_score_tbl, file = "benchmarking/summary/cell_counts_inmf.txt", sep = "\t")




########################################
############## SA-GPU #################
########################################
usages = read.csv("benchmarking/sa_gpu/output/output_W.txt", header = TRUE, sep = "\t", row.names = 1)
rownames(usages) = rownames(metaCells) 
ranks <- data.frame(gep = apply(usages, 1, which.max), score = apply(usages, 1, max))

high_score_tbl = matrix(nrow=50,ncol=40)
for (k in 1:50){
  tbl = table(factor(unlist(metaCells[rownames(metaCells) %in% rownames(usages[usages[,k]>mean(ranks[ranks$gep==k,]$score),]),1]), 
                     levels = 1:40)) 
  high_score_tbl[k,] = tbl
}
write.table(high_score_tbl, file = "benchmarking/summary/cell_counts_sagpu.txt", sep = "\t")




########################################
############## ScVI #################
########################################
usages = read.table("benchmarking/scvi/output/scvi_representation_42.txt", sep = " ", header = 0)
rownames(usages) = rownames(metaCells) 
ranks <- data.frame(gep = apply(usages, 1, which.max), score = apply(usages, 1, max))

high_score_tbl = matrix(nrow=42,ncol=40)
for (k in 1:42){
  tbl = table(factor(unlist(metaCells[rownames(metaCells) %in% rownames(usages[usages[,k]>mean(ranks[ranks$gep==k,]$score),]),1]), 
                     levels = 1:40)) 
  high_score_tbl[k,] = tbl
}
write.table(high_score_tbl, file = "benchmarking/summary/cell_counts_scvi.txt", sep = "\t")
