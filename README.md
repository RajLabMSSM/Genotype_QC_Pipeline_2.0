# Genotype QC Pipeline 2.0

Author: Ben Muller

## How to Run
### Step 1:
Move or symlink all vcf files into one folder. This is your input folder

### Step 2:
Edit `config.yaml` as necessary: 
- inFolder - the input folder containing all input files.
- suffix - options are 'vcf' or 'plink' depending on the input file type. Keep as 'vcf' for now, support for plink input will come soon!
- outFolder - the output folder where you want the output to go.
- data_name - name your dataset. All files will be concatted into one with this name. The final output will consist of files named <data_name>.MAF.vcf.gz and the like.

All other configurations allow you to keep intermediate files and change filtering thresholds for QC. They are all defined inside `config.yaml`

### Step 3:
Edit `cluster.yaml` file as necessary. This file is used to specfiy resources when submitting a snakemake job to the cluster.
Tip: Really only worry about `queue`, `mem`, and `time`. You could probably leave the whole thing as is.

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
