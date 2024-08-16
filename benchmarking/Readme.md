# Walkthrough with example data 
**Dependencies of each method is provided in the seperate requirements.txt files**.
## Prepare count data with only highly variable genes
   ```
   python benchmarking/subset_hvg.py --counts "simulated_data/counts.npz" --num_highvar_genes 2000
   ```
   Resulting file can be found on [Open Science Framework](https://osf.io/tknm2/).
   
## Run CSI-GEP from start to finish
   ```
   bsub < CSI-GEP/csigep_submit.bsub # rescue step is disabled for the complexity of the dataset
   ```
   CSI-GEP automatically determines 41 GEPs at JL = 20 & k = 47.\
   Final results can be found at CSI-GEP/output_for_benchmarking/GEP_results/.

## Run SA-GPU
   ```
   python benchmarking/sa_gpu/cmd_SA_GPU.py
   ```
   SA-GPU automatically picks k = 50. logfile can be found at benchmarking/sa_gpu/output/.\
   Output files can be found on [Open Science Framework](https://osf.io/tknm2/). 

## Run ScVI
   ```
   python benchmarking/scvi/cmd_scvi.py # comment out lines 17 & 18, uncomment line 12
   ```
   Based on the output log file, marginal ll = -185.0074 at k = 42 is relatively the lowest among reasonably large value of k
   ```
   python benchmarking/scvi/cmd_scvi.py # uncomment lines 17 & 18, comment out line 12
   ```
   Output files can be found on [Open Science Framework](https://osf.io/tknm2/) and benchmarking/scvi/output/.

## Run iNMF
   ```
   bsub < benchmarking/inmf/suggestK.bsub
   Rscript becnmarking/inmf/elbow.R # k = 52
   ```
   Modify k value in benchmarking/inmf/inmf_run.R
   ```
   bsub < benchmarking/inmf/inmf.bsub
   ```
   Output files can be found on [Open Science Framework](https://osf.io/tknm2/) and benchmarking/inmf/output/. 
   
## Summary
   1. ### Heatmaps 
      The cell usage scores are converted into contingency table by
      1. first assigning a GEP ID to each cell for which it has the highest usage score; and
      2. selecting cells only expressing a GEP with usage score higher than the avergae of all cells assigned to this GEP.
      Details can be found in the Methods section of our paper. \
      The contingency tables are reordered so that the maximum cell counts lie on the diagonal. \
      [Alt text](https://www.github.com/geeleherlab/CSI-GEP/main/benchmarking/heatmaps.pdf)

   3. ### Customized metrics
   4. ### Simplified silhouette score
   5. ### Mutual information
   6. ### Adjusted Rand index
