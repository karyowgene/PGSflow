configfile: "./config/config.yaml"

# Define samples dynamically from raw FASTQ files
samples = glob_wildcards(f"{config['raw_fastq_path']}/{{sample}}_1.fastq.gz").sample

rule all:
    input:
        expand("./output/prs_scores/{sample}_prs_report.pdf", sample=samples),
        expand("./output/single_snp_risk/{sample}_single_snp.txt", sample=samples)

rule fastqc:
    input:
        r1 = f"{config['raw_fastq_path']}/{{sample}}_1.fastq.gz",
        r2 = f"{config['raw_fastq_path']}/{{sample}}_2.fastq.gz"
    output:
        html = directory("./output/fastqc/{sample}_fastqc"),
        zip = "./output/fastqc/{sample}_fastqc.zip"
    log: "./logs/fastqc/{sample}.log"
    threads: 2
    shell:
        "fastqc {input.r1} {input.r2} -o {output.html} -t {threads} > {log} 2>&1"

rule trim:
    input:
        r1 = f"{config['raw_fastq_path']}/{{sample}}_1.fastq.gz",
        r2 = f"{config['raw_fastq_path']}/{{sample}}_2.fastq.gz"
    output:
        trim_r1 = "./output/trimmed/{sample}_1_trimmed.fastq.gz",
        trim_r2 = "./output/trimmed/{sample}_2_trimmed.fastq.gz",
        unpaired = "./output/trimmed/{sample}_unpaired.fastq.gz"
    log: "./logs/trim/{sample}.log"
    params:
        adapters = config["trimmomatic_adapters"]
    threads: 4
    shell:
        "trimmomatic PE -threads {threads} {input.r1} {input.r2} "
        "{output.trim_r1} {output.unpaired} "
        "{output.trim_r2} {output.unpaired} "
        "ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 "
        ">> {log} 2>&1"

rule align:
    input:
        r1 = "./output/trimmed/{sample}_1_trimmed.fastq.gz",
        r2 = "./output/trimmed/{sample}_2_trimmed.fastq.gz"
    output:
        bam = "./output/aligned/{sample}.bam"
    log: "./logs/align/{sample}.log"
    params:
        ref = config["reference_genome"]
    threads: 8
    shell:
        "bwa mem -t {threads} {params.ref} {input.r1} {input.r2} | "
        "samtools sort -@ {threads} -o {output.bam} - > {log} 2>&1"

rule call_variants:
    input:
        bam = "./output/aligned/{sample}.bam"
    output:
        vcf = "./output/variants/{sample}.vcf"
    log: "./logs/variants/{sample}.log"
    params:
        ref = config["reference_genome"]
    shell:
        "bcftools mpileup -f {params.ref} {input.bam} | "
        "bcftools call -mv -Ov -o {output.vcf} > {log} 2>&1"

rule vcf_to_plink:
    input:
        vcf = "./output/variants/{sample}.vcf"
    output:
        plink = "./data/user_geno/{sample}.bed"
    shell:
        "plink --vcf {input.vcf} --make-bed --out ./data/user_geno/{wildcards.sample}"

rule download_gwas_snps:
    output:
        snps = "./data/gwas_bmd_snps.tsv"
    params:
        trait = "bone mineral density",
        pval = "5e-8"
    script:
        "./scripts/query_gwas_api.py"

rule download_pgs_weights:
    output:
        weights = "./data/pgs_bmd_weights.txt"
    params:
        pgs_id = "PGS000013"
    script:
        "./scripts/query_pgs_catalog.py"

rule qc_imputation:
    input:
        plink = "./data/user_geno/{sample}.bed",
    output:
        qc_geno = "./output/qc/{sample}_qc.bed",
    shell:
        "plink2 --bfile {input.plink} --maf 0.01 --hwe 1e-6 --mind 0.1 --make-bed --out {output.qc_geno}"

rule calculate_prs:
    input:
        geno = "./output/qc/{sample}_qc.bed",
        weights = rules.download_pgs_weights.output.weights
    output:
        prs_raw = "./output/prs_scores/{sample}_raw_prs.txt"
    wildcard_constraints:
        sample=r"\w+"
    shell:
        "plink --bfile {input.geno} --score {input.weights} 1 3 2 header --out {output.prs_raw}"

rule normalize_prs:
    input:
        prs_raw = "./output/prs_scores/{sample}_raw_prs.txt",
        ref_panel = config["reference_panel"]
    output:
        prs_zscore = "./output/prs_scores/{sample}_prs_zscore.txt"
    script:
        "./scripts/prs_calculation.R"

rule single_snp_risk:
    input:
        geno = "./output/qc/{sample}_qc.bed",
        snp_list = rules.download_gwas_snps.output.snps
    output:
        report = "./output/single_snp_risk/{sample}_single_snp.txt"
    shell:
        "python ./scripts/risk_classifier.py --geno {input.geno} --snps {input.snp_list} --output {output.report}"

rule generate_report:
    input:
        prs_zscore = "./output/prs_scores/{sample}_prs_zscore.txt",
        single_snp = "./output/single_snp_risk/{sample}_single_snp.txt"
    output:
        report = "./output/prs_scores/{sample}_prs_report.pdf"
    shell:
        "Rscript ./scripts/generate_report.R {input.prs_zscore} {input.single_snp} {output.report}"
