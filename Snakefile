import pandas as pd
from snakemake.utils import validate


##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")
SAMPLES = samples.index.values.tolist()


##### rules #####

rule all:
    input:
        expand("results/{sample}/sleuth/dataframe.tsv", sample=SAMPLES)

rule index:
    input:
        config["ref_transcriptome"]
    output:
        "temp/transcripts.idx"
    shell:
        "kallisto index -i {output} {input}"

rule counts:
    input:
        tra = "temp/transcripts.idx",
        fq1 = lambda wildcards: samples['fq1'][wildcards.sample],
        fq2 = lambda wildcards: samples['fq2'][wildcards.sample]
    output:
        directory("results/{sample}/kallisto")
    conda:
        "envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.tra} -o {output} -b 10 {input.fq1} {input.fq2}"


rule sleuth:
    input:
        kal = expand("results/{sample}/kallisto", sample=SAMPLES),
        sam = config["samples"]
    output:
        expand("results/{sample}/sleuth/dataframe.tsv", sample=SAMPLES)
    conda:
        "envs/sleuth.yaml"
    script:
        "scripts/sleuth.R"

#rule boxenplot:
