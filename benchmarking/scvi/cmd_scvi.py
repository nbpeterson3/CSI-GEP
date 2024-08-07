import os
home_dir = "benchmarking/scvi"
output_dir = home_dir + "/output"
os.system('mkdir ' + output_dir)
ngenes = 2000

#### run SignatureAnalyzer
n_layer = 2

filtered_data = "benchmarking/filtered_counts_%d.npz" % (ngenes)
# search for K
scvi_cmd = "bsub -P XL -J scvi -q gpu -gpu \"num=1/host\" -R \"rusage[mem=500GB] \" -oo %s -eo %s 'python benchmarking/scvi/trainer.py --counts %s --layers %d --option 1'" % (output_dir, output_dir, filtered_data, n_layer)

# pick K based on resulting marginal log-likelihood values 

# fixed K
# K0 = 25
# scvi_cmd = "bsub -P XL -J scvi -q gpu -gpu \"num=1/host\" -R \"rusage[mem=500GB]\" -oo %s -eo %s 'python benchmarking/scvi/trainer.py --counts %s --K0 %d --layers %d --option 2'" % (output_dir, output_dir, filtered_data, K0, n_layer)


os.system(scvi_cmd)
