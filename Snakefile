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
        config['graph_p_values'],
        config['graph_counts'],
        config['graph_pca'],
        config['graph_stripplot_lowestpval']


rule index:
    input:
        config["ref_transcriptome"]
    output:
        config["transcripts_index"] # create .idx file
    conda:
        "envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input}"

rule counts:
    input:
        tra = config["transcripts_index"],
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
        #sleuth = expand("results/{sample}/sleuth/dataframe.tsv", sample=SAMPLES),
        counts_norm = 'temp/counts_normalized.tsv',
        sleuth_res = 'temp/sleuth_table.tsv',
        sleuth_matrix = 'temp/sleuth_matrix.tsv',
        sleuth_object = 'temp/sleuth_object'
    conda:
        "envs/sleuth.yaml"
    script:
        "scripts/sleuth.R"

rule boxenplot:
    input:
        sam = config["samples"],
        sle = 'temp/counts_normalized.tsv'
    output:
        config['graph_counts']
    conda:
        "envs/python_plots.yaml"
    script:
        "scripts/boxenplot.py"

rule p_values:
    input:
        'temp/sleuth_table.tsv'
    output:
        config['graph_p_values']
    conda:
        "envs/python_plots.yaml"
    script:
        "scripts/p_values.py"

rule stripplot:
    input:
        'temp/sleuth_table.tsv',
        'temp/counts_normalized.tsv'
    output:
        config['graph_stripplot_lowestpval']
    conda:
        "envs/python_plots.yaml"
    script:
        "scripts/stripplot_top20.py"

rule pca:
    input:
        sleuth_object = 'temp/sleuth_object'
    output:
        graph_pca = config['graph_pca']
    conda:
        "envs/sleuth.yaml"
    script:
        "scripts/pca.R"
