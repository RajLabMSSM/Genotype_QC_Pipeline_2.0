# Genotype QC Pipeline

Author: Ben Muller

## How to Run
### Step 1:
Move or symlink all vcf files into one folder. This is your input folder

### Step 2:
Edit `config.yaml` as necessary: 
- inFolder - input folder containing all vcf files
- suffix - keep as .vcf.gz for now
- outFolder - output folder where you want the cleaned data to be stored
- data_name - name your dataset. All files will be concatted into one with this name.

All other configurations allow you to keep intermediate files and change filtering thresholds for QC. They are all defined inside `config.yaml`

### Step 3:
Edit `cluster.yaml` file as necessary. This file is used to specfiy resources when submitting a snakemake job to the cluster.
Tip: Really only worry about `queue`, `mem`, and `time`

### Step 4:
Running the Pipeline

Always do a dry run first:
```
conda activate snakemake
snakemake -s Snakefile -npr
```
### Run snakemake in one of 3 ways:

You can run in an interactive job:

1. Like so:
```
snakemake -s Snakefile
```
 Or you can submit Snakejobs to the cluster (see `RajLabMSSM/RNA-pipelines/snakejob/` for more info).

2. Runs locally but submits jobs to cluster
```
snakejob -s Snakefile -c config.yaml
```
3. Runs and submits job on cluster
```
snakejob_HPC -s Snakefile -c config.yaml
```

## Important Notes
- Code assumes `nochr` chromosome naming.
- Lifting step not included. Must lift to hg38 before running
- input must be `.vcf` format.

Will add fixes for these (and more) later on.

## QC Pipeline Output
The pipeline will create a folder named after the `data_name` configuration specified by the user. Inside will be three folders:

1. `filtered_data` - contains two QCed data files, one including the MAF step (.MAF) and one not including the MAF step (.noMAF).
2. `relatednaess` - contains an html file displaying relatedness for manual filtering by user.
3. `population_struct` - contains an html file displaying ancestry for manual filtering by user.
