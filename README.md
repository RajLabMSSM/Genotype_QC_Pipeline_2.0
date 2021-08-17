# Genotype QC Pipeline 2.0

Author: Ben Muller

## How to Run
### Step 1:
Move or symlink all input files into one folder. This is the input folder you specify in your config.yaml file.

Input files should be in vcf format, one file per chromosome support for one big VCF will be added later. Pipeline requires build hg38 and chr prefix.

### Step 2:
Edit the `config.yaml` file as necessary: 
- in_folder - the input folder containing all input files as specified in Step 1.
- input_type - options are 'vcf' or 'plink' depending on the input file type. Keep as 'vcf' for now, support for plink input will come soon!
- out_folder - the output folder where you want the output to go.
- data_name - name your dataset. All chrom files will be concatted into one big file with this name. The final output will consist of files named <data_name>.MAF.vcf.gz and the like.

All other configurations allow you to change filtering thresholds for QC and keep intermediate files for debugging. They are all defined inside `config.yaml`.

### Step 3:
Edit the `cluster.yaml` file as necessary. This file is used to specify minerva resources when submitting a snakemake job to the cluster.
Tip: Really only worry about `queue`, `mem`, and `time`. You could probably leave the whole thing as is.

### Step 4:
Running the Pipeline.

Always do a dry run first:
```
conda activate snakemake
snakemake -s Snakefile -npr
```
### Run snakemake in one of 3 ways:

You can run in an interactive job:

1. Like so:
```
snakemake -s Snakefile --cores <num_cores>
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
- As mentioned above, the pipeline assumes `chr` chromosome naming.
- A liftover step not included. you must lift to hg38 before running
- input must be `.vcf.gz` format, with one file per chrom.

I'll add fixes for these (and more) later on.

## QC Pipeline Output`
There will be 3 main output folders:

1. `clean_data` - contains two QCed data files, one including the MAF step (.MAF) and one not including the MAF step (.noMAF).
2. `population` - contains an html and tsv files containing information on both relatedness (`somalier.*`) and ancestry (`somalier-ancestry.somalier-ancestry.*`). Use these to manually filter your data as you see fit.
3. `QC_stats` - contains filtering statistics for each of the 3 QC steps: missingness, MAF, and Hardy Weinberg.



