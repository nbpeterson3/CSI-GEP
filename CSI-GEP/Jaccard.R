args = (commandArgs(trailingOnly = TRUE))

library(jaccard)
library(igraph)
library(vroom)
library(foreach)
library(doParallel)
registerDoParallel(16)

set.seed(1234)

#This function computes the Jaccard index and significance test for all programs in each benchmark split
#Once this is completed, a network graph is built using igraph package
#Community detection of the network graph is then calculated using betweeness scores

calcSplitStats <- function(split0, split1, k, jval){
  print("calcSplitStats")
  # for each row of split0, compare it to every row of split1
  d <-
  foreach(i = 1:ncol(split0), .combine = 'rbind') %:%
    foreach(j = 1:ncol(split1), .combine = 'rbind') %dopar% {
      GEP_split0 = as.numeric(rank(-split0[, i]) < jval)
      GEP_split1 = as.numeric(rank(-split1[, j]) < jval)
      jaccardTest <- jaccard.test.mca(GEP_split0, GEP_split1)
      data.frame(Ind = jaccardTest$statistics, IndPs = jaccardTest$pvalue, split0 = i, split1 = j)
    }
  fwerJaccard <- p.adjust(d$IndPs, method="bonferroni")
  
  # Calculate FWERs
  k <- as.numeric(k) 
  # fwerJaccard <- p.adjust(jaccardIndPs, method="bonferroni")
  fwerJaccard <- t(matrix(fwerJaccard, nrow = k, ncol = k))
  
  # Calculate number of significant findings.
  sigFwerJaccard <- which(fwerJaccard < 0.05, arr.ind = T)
  
  #Calculate community number
  if(nrow(sigFwerJaccard) > 0){
    community.df <- as.data.frame(sigFwerJaccard)
    community.df$row <- paste0("Split0_GEP", community.df$row)
    community.df$col <- paste0("Split1_GEP", community.df$col)
    community.graph <- graph_from_data_frame(community.df, directed = F)
    ceb <- cluster_edge_betweenness(community.graph)
    community_number <- length(communities(ceb))
  }else{
    community_number <- 0
  }
  return(list(jaccardIndPs=d$IndPs, jaccardInd=d$Ind, sigFwerJaccard=sigFwerJaccard, jaccard.com=community_number))
}

#Obtain the data directory and rank for which the calcSplitStats function should be run              
k <- args[1]                                                                       
data.dir <- args[2]

#jaccard lengths
jaccard.values <- c(seq(10,100,10))

#Results list
#Gets written to RDS file at end of R script

#Read in NMF gene spectra scores  
print(paste0("Rank: ", k))
split0 <- t(vroom(file = paste0(data.dir, "cNMF_split0/cNMF_split0.gene_spectra_score.k_", k, ".dt_0_1.txt"), col_select = c(-...1), show_col_types = FALSE))
split1 <- t(vroom(file = paste0(data.dir, "cNMF_split1/cNMF_split1.gene_spectra_score.k_", k, ".dt_0_1.txt"), col_select = c(-...1), show_col_types = FALSE))
print(paste0("split0 Filename: cNMF_split0.gene_spectra_score.k_", k, ".dt_0_1.txt"))
print(paste0("split1 Filename: cNMF_split1.gene_spectra_score.k_", k, ".dt_0_1.txt"))

#Compute jaccard index for different jaccard lengths
foreach (jval = 1:length(jaccard.values)) %dopar% {
  print(paste0("Jaccard Length: ", jaccard.values[jval]))
  splitstatlist <- list()
  splitstatlist <- calcSplitStats(split0, split1, k, jaccard.values[jval]) 
  jval.name <- as.character(jaccard.values[jval])
  saveRDS(splitstatlist, file = paste0(data.dir, "results/jaccardtest_results_K", k, "_JL", jval.name, ".RDS"))
}

jaccardtest_results <- list()

for (j in jaccard.values){
  jval.name <- as.character(j)
  splitstatlist = readRDS(paste0(data.dir, "results/jaccardtest_results_K", k, "_JL", jval.name, ".RDS"))
  jaccardtest_results[[jval.name]] <- splitstatlist
}
#This list contains the computations across all jaccard lengths for the rank provided on the command line
saveRDS(jaccardtest_results, file = paste0(data.dir, "results/jaccardtest_results_K", k, ".RDS"))

