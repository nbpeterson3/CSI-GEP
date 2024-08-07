import os
home_dir = "benchmarking/sa_gpu"
output_dir = home_dir + "/output"
os.system('mkdir ' + output_dir)
K_init = 100
ngenes = 2000


#### run SignatureAnalyzer
filtered_data = "benchmarking/filtered_counts_%d.npz" % (ngenes)

SA_cmd = "bsub -P XL -J SA_GPU -q gpu -gpu \"num=1/host\" -R \"rusage[mem=500GB]\" -oo %s -eo %s 'python benchmarking/sa_gpu/SignatureAnalyzer-GPU.py --data %s --max_iter=100000 --output_dir %s --K0 %d --prior_on_W L1 --prior_on_H L2 --labeled'" % (output_dir, output_dir, filtered_data, output_dir, K_init)
os.system(SA_cmd)
