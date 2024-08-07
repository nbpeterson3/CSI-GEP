library(rliger)
library(Matrix)
library(data.table)
library(reticulate)

np = import("numpy")
pd = import("pandas")
data = np$load("benchmarking/filtered_counts_2000.npz", allow_pickle = TRUE)
counts = data$f[["data"]]
input = pd$DataFrame(counts, index = data$f[["index"]], columns = data$f[["columns"]])


obj = createLiger(list(HTAPP = t(input)))
obj = normalize(obj)
obj = selectGenes(obj)
obj = scaleNotCenter(obj, remove.missing = F)
obj = online_iNMF(obj, k = 60, miniBatch_size = 5000, max.epochs = 5, verbose = FALSE)

write.table(obj@W, 
    paste0("benchmarking/inmf/output/iNMF_spectra.k_60.txt"))

write.table(obj@H, 
    paste0("benchmarking/inmf/output/iNMF_usages.k_60.txt"))

