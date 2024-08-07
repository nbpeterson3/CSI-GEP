# CSI-GEP
Consensus and Scalable Interpretations of Gene Expression Programs in Single Cells

## Benchmarking 
## Docker
We provide Docker images that minimizes conplexity of configuration, and a step-by-step general instruction on running it on HPC.
1. On your local host, pull Docker images:
   ```
   docekr pull ghcr.io/geeleherlab/csi-gep_py:latest
   docekr pull ghcr.io/geeleherlab/csi-gep_r:latest
   ```
2. Push to the HPC system of your organization. Below is an example on St. Jude HPC that utilizes singularity:
   ```
   singularity remote login --user <my_userID> docker://svlprhpcreg01.stjude.org
   docker image push svlprhpcreg01.stjude.org/hpcf/<docker_image_name>
   ```
3. Submit jobs. An example code csigep_submit.bsub can be found at Docker/tutorial/.
