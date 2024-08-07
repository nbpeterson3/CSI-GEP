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
   CSI-GEP automatically determines 41 GEPs at JL = 20 & k = 47

## Run SA-GPU
   ```
   python benchmarking/sa_gpu/cmd_SA_GPU.py
   ```
   SA-GPU automatically picks k = 50

## Run ScVI
   ```
   python benchmarking/scvi/cmd_scvi.py # comment out lines 17 & 18, uncomment line 12
   ```
   Based on the output log file, marginal ll = 245.57 at k = 25 is relatively the lowest among reasonably large value of k
   ```
   python benchmarking/scvi/cmd_scvi.py # uncomment lines 17 & 18, comment out line 12
   ```

## Run iNMF
   ```
   bsub < benchmarking/inmf/suggestK.bsub
   Rscript becnmarking/inmf/elbow.R # k = 60
   ```
   Modify k value in benchmarking/inmf/inmf_run.R
   ```
   bsub < benchmarking/inmf/inmf.bsub
   ```
## Summary
   1. ### Heatmaps
   2. ### Customized metrics
   3. ### Simplified silhouette score
   4. ### Mutual information
   5. ### Adjusted Rand index
