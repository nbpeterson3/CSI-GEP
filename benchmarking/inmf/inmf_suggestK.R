library(rliger)
library(foreach)
library(doParallel)
library(Matrix)
library(data.table)
#registerDoParallel(parallel::detectCores()/2)
# registerDoParallel(16)
library(reticulate)
library(ggplot2)

dir.create(file.path("benchmarking/inmf/", "output"), showWarnings = FALSE)

np = import("numpy")
pd = import("pandas")
data = np$load("benchmarking/filtered_counts_2000.npz", allow_pickle = TRUE)
counts = data$f[["data"]]
input = pd$DataFrame(counts, index = data$f[["index"]], columns = data$f[["columns"]])


obj = createLiger(list(HTAPP = t(input))) 
obj = normalize(obj)
obj = selectGenes(obj)
obj = scaleNotCenter(obj)


kl_divergence_uniform = function(object, Hs=NULL)
{
  if (is.null(Hs)) {Hs = object@H}
  n_cells = sum(sapply(Hs, nrow))
  n_factors = ncol(Hs[[1]])
  dataset_list = list()
  for (i in 1:length(Hs)) {
    scaled = scale(Hs[[i]], center=FALSE, scale=TRUE)
    
    inflated = t(apply(scaled, 1, function(x) {
      replace(x, x == 0, 1e-20)
    }))
    inflated = inflated/rowSums(inflated)
    divs = apply(inflated, 1, function(x) {log2(n_factors) + sum(log2(x) * x)})
    dataset_list[[i]] = divs
  }
  return(dataset_list)
}

k.test = c(seq(2,60), seq(70,200,10))
verbose = TRUE
nrep = 1
rand.seed = 1
lambda = 5
thresh = 1e-4
max.iters = 100
gen.new = FALSE
num.cores = 16

time_start <- Sys.time()
# optimize largest k value first to take advantage of efficient updating
if (verbose) {
  message("This operation may take several minutes depending on number of values being tested")
}
rep_data <- list()
for (r in 1:nrep) {
  if (verbose) {
    message("Preprocessing for rep ", r,
            ": optimizing initial factorization with largest test k=",
            k.test[length(k.test)])
  }
  object <- optimizeALS(obj, k = k.test[length(k.test)], lambda = lambda, thresh = thresh,
                        max.iters = max.iters, nrep = 1, rand.seed = (rand.seed + r - 1))
  if (verbose) {
    message('Testing different choices of k')
  }
  cl <- parallel::makeCluster(num.cores)
  doParallel::registerDoParallel(cl)
  
  i <- 0
  data_matrix <- foreach(i = length(k.test):1, .combine = "rbind",
                         .packages = 'rliger') %dopar% {
                           if (i != length(k.test)) {
                             if (gen.new) {
                               ob.test <- optimizeALS(object,
                                                      k = k.test[i], lambda = lambda, thresh = thresh,
                                                      max.iters = max.iters, rand.seed = (rand.seed + r - 1)
                               )
                             } else {
                               ob.test <- optimizeNewK(object,
                                                       k.new = k.test[i], lambda = lambda, thresh = thresh,
                                                       max.iters = max.iters, rand.seed = (rand.seed + r - 1)
                               )
                             }
                           } else {
                             ob.test <- object
                           }
                           dataset_split <- kl_divergence_uniform(ob.test)
                           unlist(dataset_split)
                         }
  #close(pb)
  parallel::stopCluster(cl)
  data_matrix <- data_matrix[nrow(data_matrix):1, ]
  rep_data[[r]] <- data_matrix
}

medians <- Reduce(cbind, lapply(rep_data, function(x) {apply(x, 1, median)}))
if (is.null(dim(medians))) {
  medians <- matrix(medians, ncol = 1)
}
mean_kls <- apply(medians, 1, mean)

time_elapsed <- difftime(Sys.time(), time_start, units = "auto")
if (verbose) {
  cat(paste("\nCompleted in:", as.double(time_elapsed), units(time_elapsed)))
}
# make dataframe
df_kl <- data.frame(median_kl = c(mean_kls, log2(k.test)), k = c(k.test, k.test),
                    calc = c(rep('KL_div', length(k.test)), rep('log2(k)', length(k.test))))

print(df_kl)
write.csv(df_kl[df_kl$calc == "KL_div",], "benchmarking/inmf/output/median_kl.csv")

p1 <- ggplot(df_kl, aes_string(x = 'k', y = 'median_kl', col = 'calc')) + geom_line(size=1) +
    geom_point() +
    theme_classic() + labs(y='Median KL divergence (across all cells)', x = 'K') +
    guides(col=guide_legend(title="", override.aes = list(size = 2))) +
    theme(legend.position = 'top')

pdf("benchmarking/inmf/output/suggestK_plot.pdf")
print(p1)
dev.off()



