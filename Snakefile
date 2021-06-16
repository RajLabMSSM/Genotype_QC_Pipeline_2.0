import os
import glob
import pathlib

configfile: "config.yaml"

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

# Make intermediate directories here
pathlib.Path(out_folder + "prepped_data").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "missing").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "MAF").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "hwe").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "population").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "clean_data").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "QC_stats").mkdir(parents=True, exist_ok=True)
pathlib.Path(out_folder + "temp").mkdir(parents=True, exist_ok=True)

rule all:
    input:
        #temp_output
        expand(out_folder + 'MAF/{file_name}.vcf.gz', file_name=file_names),
        expand(out_folder + "clean_data/" + data_name + "{suff}", suff=['_MAF.vcf.gz', '_noMAF.vcf.gz']),
        #out_folder + "clean_data/" + data_name + "_noMAF.vcf.gz",
        #out_folder + "population/somalier.html",
        out_folder + "population/somalier-ancestry.somalier-ancestry.html"
        #ancestry_output,
        #relatedness_output,
        #outFolder + data_name + "/filtered_data/" + data_name + ".MAF.vcf.gz",
        #outFolder + data_name + "/filtered_data/" + data_name + ".noMAF.vcf.gz"

rule prep_data:
    input: 
        in_folder + "{file_name}" + suffix
    output: 
        expand(out_folder + "prepped_data/{{file_name}}{plink_suffix}", plink_suffix=plink_suffix)
    run:
        if file_type == "vcf":
            shell("ml plink/2.0dev060421;\
                   plink2 --vcf {input} --max-alleles 2 --make-pgen --out {out_folder}prepped_data/{wildcards.file_name};")

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
               plink2 --pfile {out_folder}hwe/{wildcards.file_name}_step2 --recode vcf "bgz" --out {out_folder}hwe/{wildcards.file_name}')

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
               plink2 --pfile {out_folder}hwe/{wildcards.file_name}_step2 --maf {params.MAF} --recode vcf "bgz" --out {out_folder}MAF/{wildcards.file_name}')

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

rule filter_files:
    input:
        out_folder + "temp/sample_intersection.txt"
        #expand(out_folder + "temp/{file_name}_samples.txt", file_name=file_names)
        #out_folder + "hwe/{file_name}.vcf.gz",
        #out_folder + "MAF/{file_name}.vcf.gz"
    output:
        out_folder + "temp/{file_name}_MAF.vcf.gz",
        out_folder + "temp/{file_name}_noMAF.vcf.gz"
    run:
        #shell('ml bcftools;\
        #       rm -rf {out_folder}temp; mkdir -p {out_folder}temp;\
        #       bcftools query -l {out_folder}/MAF/{wildcards.file_name}.vcf.gz >> {out_folder}/temp/{wildcards.file_name}_samples.txt;\
        shell('ml bcftools;\
               bcftools view --force-samples -S {input} {out_folder}MAF/{wildcards.file_name}.vcf.gz -Oz > {out_folder}/temp/{wildcards.file_name}_MAF.vcf.gz;\
               bcftools view --force-samples -S {input} {out_folder}hwe/{wildcards.file_name}.vcf.gz -Oz > {out_folder}/temp/{wildcards.file_name}_noMAF.vcf.gz')


rule merge_files:
    input:
        expand(out_folder + "temp/{file_name}_MAF.vcf.gz", file_name=file_names),
        expand(out_folder + "temp/{file_name}_noMAF.vcf.gz", file_name=file_names)
    output:
        out_folder + "clean_data/" + data_name + "{suff}"
        #out_folder + "clean_data/" + data_name + "_noMAF.vcf.gz",
        #out_folder + "clean_data/" + data_name + "_MAF.vcf.gz"
    run:
        #shell('ml bcftools;\
               #rm -rf {out_folder}temp1; mkdir -p {out_folder}temp1;\
               #for i in {out_folder}/MAF/*gz; do bcftools query -l $i >> {out_folder}/temp1/$(basename -s .vcf.gz $i)_samples.txt; done;\
               #sort {out_folder}/temp1/* | uniq -d >> {out_folder}/temp1/sample_intersection.txt;\
               #for i in {file_names}; do bcftools view --force-samples -S {out_folder}/temp1/sample_intersection.txt {out_folder}/MAF/$i.vcf.gz -Oz > {out_folder}/temp1/$i.MAF.vcf.gz; bcftools view --force-samples -S {out_folder}temp1/sample_intersection.txt {out_folder}hwe/$i.vcf.gz -Oz > {out_folder}/temp1/$i.noMAF.vcf.gz; done;\
        shell('ml bcftools;\
               bcftools concat {out_folder}/temp/*{wildcards.suff} | bgzip -c > {out_folder}clean_data/{data_name}{wildcards.suff};\
               tabix {out_folder}clean_data/{data_name}{wildcards.suff}')
               #bcftools concat {out_folder}/temp/*.MAF.vcf.gz | bgzip -c > {out_folder}clean_data/{data_name}_MAF.vcf.gz;\
               #bcftools concat {out_folder}/temp/*.noMAF.vcf.gz | bgzip -c > {out_folder}clean_data/{data_name}_noMAF.vcf.gz;\
               #tabix {out_folder}clean_data/{data_name}_noMAF.vcf.gz')


#rule merge_convert_to_vcf:
#    input:
#        expand(out_folder + "hwe/{file_name}_step2{plink_suffix}", file_name=file_names, plink_suffix=plink_suffix),
#        expand(out_folder + "MAF/{file_name}{plink_suffix}", file_name=file_names, plink_suffix=plink_suffix)
#    output:
#        out_folder + "clean_data/" + data_name + "_MAF.vcf.gz",
#        out_folder + "clean_data/" + data_name + "_noMAF.vcf.gz"
#    run:
#        #shell('for i in {out_folder}/MAF/*; do echo $(pw)${i%.*} >> mergelist.txt; done;\
#        #       plink2 --pmerge-list mergelist.txt --sort-vars --make-pgen --out')
#        shell('ml plink/2.0dev060421;\
#               mkdir {out_folder}/clean_data/temp1; mkdir {out_folder}/clean_data/temp2;\
#               for i in {file_names};\
#               do plink2 --pfile {out_folder}MAF/$i --recode vcf "bgz" --out {out_folder}/clean_data/temp1/$i_MAF\
#               do plink2 --pfile {out_folder}hwe/$i_step2 --recode vcf "bgz" --out {out_folder}/clean_data/temp1/$i_noMAF'\
#               )

rule relatedness_ancestry:
    input:
        out_folder + "clean_data/" + data_name + "_noMAF.vcf.gz"
    output:
        out_folder + "population/somalier.html",
        out_folder + "population/somalier-ancestry.somalier-ancestry.html"
    params:
        sites="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/sites/sites.hg38.nochr.vcf.gz",
        ref="/sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.fa",
        ancestry_labels="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/ancestry/ancestry-labels-1kg.tsv",
        labeled_samples="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/ancestry/1kg-somalier/"
    run:
        shell('somalier=/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/somalier;\
               $somalier extract -d {out_folder}population/extracted/ --sites {params.sites} -f {params.ref} {out_folder}clean_data/{data_name}_noMAF.vcf.gz;\
               cd {out_folder}population/;\
               $somalier relate extracted/*.somalier;\
               $somalier ancestry --labels {params.ancestry_labels} {params.labeled_samples}*.somalier ++ extracted/*.somalier')







