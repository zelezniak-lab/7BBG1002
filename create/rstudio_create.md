# Running R Studio server on CREATE HPC

To connect to Rstudio on CREATE, first login to HPC by using terminal window.
After logging in, copy the following Rstudio starting script to your home directory:

```bash
cp /scratch/grp/msc_appbio/ops-rstudio2.0.sh .
```
I suggest you study the content inside. The script loads modules required to run working version of R and R studio on CREATE HPC.

Next, start R studio server as batch job.

```bash
sbatch ops-rstudio2.0.sh 
```
Wait a minute or two, and follow the instructions inside SLURM/sbatch job output script by using `cat  slurm-XXXXXXXX.out` where X's your last sbatch job identifier. 

In short, you need to:

- Open another terminal window and run an SSH tunnel from HPC to your local machine, using the command you will find in SLURM/sbatch job output.
- Log in to Rstudio from your internet browser using  `http://localhost:8787` using the your K number and a generated random password you will find inside SLURM/sbatch output log. 

