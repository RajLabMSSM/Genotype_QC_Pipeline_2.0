import os
import glob
import pathlib

# DON"T HARD CODE THE NAME OF THE CONFIG FILE


in_folder = config['in_folder']
file_type = config['input_type']
out_folder = config['out_folder']
data_name = config['data_name']
keep_int= config['keep_intermediate']

#geno = config['miss_SNP']
#mind = ['miss_ind']


plink_suffix = ['.psam', '.pvar', '.pgen']
suffix = ''

if file_type == 'vcf':
    suffix = '.vcf.gz'

file_paths = glob.glob(in_folder + "*" + suffix)
file_names = [os.path.basename(f).replace(suffix, "") for f in file_paths]

# all final output files
ancestry_output = out_folder + "ancestry/somalier.somalier-ancestry.html"
relatedness_output = out_folder + "relatedness/somalier.html"

# a temp output to help me debug while creating pipeline
temp_output = expand(out_folder + '{file_name}' + suffix, file_name = file_names)

print('in_folder: ', in_folder)
print('file_type: ', file_type)
print('out_folder: ', out_folder)
print('data_name: ', data_name)
print('suffix: ', suffix)
print('file_names:', file_names)
print('keep_int', keep_int, type(keep_int))

# Make intermediate directories here
pathlib.Path(out_folder + "prepped_data").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "missing").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "MAF").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "hwe").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "population").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "clean_data").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "QC_stats").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "temp").mkdir(parents=True, exist_ok=True)

#ruleorder: prep_data > missingness > HWE > MAF > generate_sample_IDs > filter_files > merge_files > index_files > sort_files > relatedness_ancestry
#wildcard_constraints:
#    suff=['_MAF.vcf.gz', '_noMAF.vcf.gz']

rule all:
    input:
        #temp_output
        #expand(out_folder + 'MAF/{file_name}.vcf.gz', file_name=file_names),
        expand(out_folder + "clean_data/" + data_name + "_anno.{suff}", suff=['MAF.vcf.gz', 'noMAF.vcf.gz']),
        expand(out_folder + "stats/" + data_name + "_anno.{suff}_stats.txt", suff=['MAF.vcf.gz', 'noMAF.vcf.gz']),
        #out_folder + "clean_data/" + data_name + "_noMAF.vcf.gz",
        #out_folder + "population/somalier.html",
        out_folder + "population/" + data_name + ".somalier.ancestry.html",
        out_folder + "population/" + data_name + ".somalier.relatedness.html",
        out_folder + "population/" + data_name + ".somalier.ancestry.PCs.tsv"
        #out_folder + "population/somalier-ancestry.somalier-ancestry.html"
        #ancestry_output,
        #relatedness_output,
        #outFolder + data_name + "/filtered_data/" + data_name + ".MAF.vcf.gz",
        #outFolder + data_name + "/filtered_data/" + data_name + ".noMAF.vcf.gz"

# Convert input data (for now only VCFs) into plink2 files - 1 file per chromosome
rule prep_data:
    input: 
        in_folder + "{file_name}" + suffix
    output: 
        expand(out_folder + "prepped_data/{{file_name}}{plink_suffix}", plink_suffix=plink_suffix)
    run:
        if file_type == "vcf":
            shell("ml plink/2.0dev060421;\
                   plink2 --vcf {input} --max-alleles 2 --make-pgen --out {out_folder}prepped_data/{wildcards.file_name};")

# Apply Missingness filters:
# First round of filtering with a higher missingness threshold to account for serious technical errors. 
# Second round uses a more standard filtering threshold
rule missingness:
    input: 
        expand(out_folder + "prepped_data/{{file_name}}{plink_suffix}", plink_suffix=plink_suffix)
    output: 
        expand(out_folder + "missing/{{file_name}}_mind2{plink_suffix}", plink_suffix=plink_suffix)
    params:
        geno = config['miss_SNP'],
        mind = config['miss_ind']
    run:
        print("geno param = ", params.geno)
        print("mind param = ", params.mind)
        shell('ml plink/2.0dev060421;\
               plink2 --pfile {out_folder}prepped_data/{wildcards.file_name} --missing --out {out_folder}QC_stats/{wildcards.file_name}_missingness;\
               plink2 --pfile {out_folder}prepped_data/{wildcards.file_name} --geno .2 --make-pgen --out {out_folder}missing/{wildcards.file_name}_geno1;\
               plink2 --pfile {out_folder}missing/{wildcards.file_name}_geno1 --mind .2 --make-pgen --out {out_folder}missing/{wildcards.file_name}_mind1;\
               plink2 --pfile {out_folder}missing/{wildcards.file_name}_mind1 --geno {params.geno} --make-pgen --out {out_folder}missing/{wildcards.file_name}_geno2;\
               plink2 --pfile {out_folder}missing/{wildcards.file_name}_geno2 --mind {params.mind} --make-pgen --out {out_folder}missing/{wildcards.file_name}_mind2;')

# Apply Hardy Weinberg filters
rule HWE:
    input:
        expand(out_folder + "missing/{{file_name}}_mind2{plink_suffix}", plink_suffix=plink_suffix)
    output:
        out_folder + "hwe/{file_name}.vcf.gz"
    params:
        step1 = config['hwe_step1'],
        step2 = config['hwe_step2']
    run:
        shell('ml plink/2.0dev060421;\
               plink2 --pfile {out_folder}missing/{wildcards.file_name}_mind2 --hardy --out {out_folder}QC_stats/{wildcards.file_name}_hwe;\
               plink2 --pfile {out_folder}missing/{wildcards.file_name}_mind2 --hwe {params.step1} --make-pgen --out {out_folder}hwe/{wildcards.file_name}_step1;\
               plink2 --pfile {out_folder}hwe/{wildcards.file_name}_step1 --hwe {params.step2} --make-pgen --out {out_folder}hwe/{wildcards.file_name}_step2;\
               plink2 --pfile {out_folder}hwe/{wildcards.file_name}_step2 --recode vcf "bgz" --output-chr chrM --out {out_folder}hwe/{wildcards.file_name}')

# Apply Minor Allele Frequency filter
# Final data will include a MAFed dataset and an unMAFed dataset
rule MAF:
    input:
        out_folder + "hwe/{file_name}.vcf.gz"
        #expand(out_folder + "hwe/{{file_name}}_step2{plink_suffix}", plink_suffix=plink_suffix)
    output:
        out_folder + "MAF/{file_name}.vcf.gz"
    params:
        MAF = config['MAF']
    run:
        shell('ml plink/2.0dev060421;\
               plink2 --pfile {out_folder}hwe/{wildcards.file_name}_step2 --freq --out {out_folder}QC_stats/{wildcards.file_name}_MAF;\
               plink2 --pfile {out_folder}hwe/{wildcards.file_name}_step2 --maf {params.MAF} --recode vcf "bgz" --output-chr chrM --out {out_folder}MAF/{wildcards.file_name}')

# Generate a list of sample IDs from each of the chromosome files.
# Create a new list containing only the intersection of all sample IDs - i.e. only the samples that exist in ALL chrom files
rule generate_sample_IDs:
    input:
        expand(out_folder + "hwe/{file_name}.vcf.gz", file_name=file_names),
        expand(out_folder + "MAF/{file_name}.vcf.gz", file_name=file_names)
    output:
        out_folder + "temp/sample_intersection.txt"
        #expand(out_folder + "temp/{file_name}_samples.txt", file_name=file_names)
    run:
        shell('ml bcftools;\
               rm -rf {out_folder}temp; mkdir -p {out_folder}temp;\
               for i in {out_folder}MAF/*.vcf.gz; do bcftools query -l $i >> {out_folder}temp/$(basename -s .vcf.gz $i)_samples.txt; done;\
               sort {out_folder}/temp/* | uniq -d >> {out_folder}/temp/sample_intersection.txt')

# Filter each chrom file using the new list of intersected IDs
rule filter_files:
    input:
        out_folder + "temp/sample_intersection.txt"
        #expand(out_folder + "temp/{file_name}_samples.txt", file_name=file_names)
        #out_folder + "hwe/{file_name}.vcf.gz",
        #out_folder + "MAF/{file_name}.vcf.gz"
    output:
        out_folder + "temp/{file_name}_sorted_MAF.vcf.gz",
        out_folder + "temp/{file_name}_sorted_noMAF.vcf.gz"
        #shell('ml bcftools;\
        #       rm -rf {out_folder}temp; mkdir -p {out_folder}temp;\
        #       bcftools query -l {out_folder}/MAF/{wildcards.file_name}.vcf.gz >> {out_folder}/temp/{wildcards.file_name}_samples.txt;\
    shell:
        'ml bcftools;\
         bcftools view --force-samples -S {input} {out_folder}MAF/{wildcards.file_name}.vcf.gz -Oz > {out_folder}/temp/{wildcards.file_name}_MAF.vcf.gz;\
         bcftools view --force-samples -S {input} {out_folder}hwe/{wildcards.file_name}.vcf.gz -Oz > {out_folder}/temp/{wildcards.file_name}_noMAF.vcf.gz;\
         tabix {out_folder}/temp/{wildcards.file_name}_MAF.vcf.gz;\
         tabix {out_folder}/temp/{wildcards.file_name}_noMAF.vcf.gz;\
         bcftools sort {out_folder}/temp/{wildcards.file_name}_MAF.vcf.gz -Oz -o {out_folder}/temp/{wildcards.file_name}_sorted_MAF.vcf.gz;\
         bcftools sort {out_folder}/temp/{wildcards.file_name}_noMAF.vcf.gz -Oz -o {out_folder}/temp/{wildcards.file_name}_sorted_noMAF.vcf.gz'

# Concatenate individual chrom files into one.
# Then sort and index VCFs. These are the final output VCF data
rule merge_files:
    input:
        #expand(out_folder + "temp/{file_name}_MAF.vcf.gz", file_name=file_names),
        #expand(out_folder + "temp/{file_name}_noMAF.vcf.gz", file_name=file_names)
        expand(out_folder + "temp/{file_name}_sorted_{{suff}}", file_name=file_names)
    output:
        out_folder + "clean_data/" + data_name + ".{suff}"
        #out_folder + "temp2/" + "merged_" + data_name + "{suff}.tbi"
    shell: 
        'ml bcftools;\
         bcftools concat {input} | bgzip -c > {out_folder}clean_data/{data_name}.{wildcards.suff}'


# New rule to index the file. 
# Why? Because if the indexing fails as it's doing now, a new rule will prevent the need to re-concat all the files which takes bloody forever.
rule index_files:
    input:
        out_folder + "clean_data/" + data_name + ".{suff}"
    output:
        out_folder + "clean_data/" + data_name + ".{suff}.tbi"
    shell:
        'ml bcftools;\
         tabix {input}'


# annotate VCF with SNP information
rule annotate_VCF:
    input:
        vcf = out_folder + "clean_data/" + data_name + ".{suff}",
        tbi =  out_folder + "clean_data/" + data_name + ".{suff}.tbi"
    output:
        vcf = out_folder + "clean_data/" + data_name + "_anno.{suff}",
        tbi = out_folder + "clean_data/" + data_name + "_anno.{suff}.tbi"
    params:
        ensembl = "/sc/arion/projects/ad-omics/data/references/hg38_reference/ensembl/ensembl_v99_hg38.vcf.gz" 
    shell:
        "ml bcftools/1.9;"
        "bcftools annotate -Oz -o {output.vcf} -a {params.ensembl} -c ID {input.vcf};"
        "tabix {output.vcf}"

# sorting vcf file
#rule sort_files:
#    input:
#        #out_folder + "temp2/" + "merged_" + data_name + "{suff}"
#        out_folder + "temp2/" + "merged_" + data_name + "{suff}.tbi"
#    output:
#        out_folder + "clean_data/" + data_name + "{suff}"
#        #out_folder + "clean_data/" + data_name + "{suff}.tbi"
#    run:
#        shell('ml bcftools;\
#               bcftools sort {out_folder}temp2/merged_{data_name}{wildcards.suff} -Oz -o {out_folder}clean_data/{data_name}{wildcards.suff};\
#               tabix {out_folder}clean_data/{data_name}{wildcards.suff}')


# Use Somalier software to generate list of relatedness and ancestry among all samples
# This step does NOT perform filtration - user should manually filter based on output
# This step DOES delete all the intermediate files though (if specified by user)
rule relatedness_ancestry:
    input:
        vcf = out_folder + "clean_data/" + data_name + "_anno.noMAF.vcf.gz",
        tbi = out_folder + "clean_data/" + data_name + "_anno.noMAF.vcf.gz.tbi"
    output:
        out_folder + "population/" + data_name + ".somalier.ancestry.html",
        out_folder + "population/" + data_name + ".somalier.relatedness.html",
        out_folder + "population/" + data_name + ".somalier.ancestry.tsv"
    params:
        sites="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/sites/sites.hg38.vcf.gz",
        ref="/sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.fa",
        ancestry_labels="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/ancestry/ancestry-labels-1kg.tsv",
        labeled_samples="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/ancestry/1kg-somalier/"
    run:
        shell('somalier=/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/somalier;\
               $somalier extract -d {out_folder}population/extracted/ --sites {params.sites} -f {params.ref} {input.vcf};\
               cd {out_folder}population/;\
               $somalier relate extracted/*.somalier;\
               $somalier ancestry --labels {params.ancestry_labels} {params.labeled_samples}*.somalier ++ extracted/*.somalier;\
               rename somalier {data_name}.somalier somalier* ;\
               mv {data_name}.somalier.html {data_name}.somalier.relatedness.html; \
               mv {data_name}.somalier-ancestry.somalier-ancestry.tsv {data_name}.somalier.ancestry.tsv; \
               mv {data_name}.somalier-ancestry.somalier-ancestry.html {data_name}.somalier.ancestry.html; \
             ')
        if keep_int==0:
            shell('rm -rf {out_folder}prepped_data;\
                   rm -rf {out_folder}missing;\
                   rm -rf {out_folder}MAF;\
                   rm -rf {out_folder}hwe;\
                   rm -rf {out_folder}temp')


rule make_PCs_for_tensorQTL:
    input:
        out_folder + "population/" + data_name + ".somalier.ancestry.tsv"
    output:
        out_folder + "population/" + data_name + ".somalier.ancestry.PCs.tsv"
    params:
        script = "scripts/make_ancestry_PCs.R"
    shell:
        "ml R/4.0.3;"
        "Rscript {params.script} --inFile {input} --outFile {output}"


rule get_stats:
    input:
        vcf = out_folder + "clean_data/" + data_name + "_anno.{suff}",
        tbi = out_folder + "clean_data/" + data_name + "_anno.{suff}.tbi"
    output:
        out_folder + "stats/" + data_name + "_anno.{suff}_stats.txt"
    shell:
        "ml bcftools; "
        "bcftools stats {input.vcf} > {output} "

