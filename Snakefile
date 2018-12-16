import pandas as pd
from snakemake.utils import validate


##### load config and sample sheets #####

configfile: "config.yaml"
# validate(config, schema="schemas/config.schema.yaml")

table = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(samples, schema="schemas/samples.schema.yaml")
SAMPLES = table.index.values.tolist()


##### rules #####

rule all:
    input:
        expand("output/abundances/{sample}", sample=SAMPLES)

rule index:
    input:
        config["ref_transcriptome"]
        #"data/ref/transcriptome.chr21.fa"
    output:
        "output/transcripts.idx"
    shell:
        "kallisto index -i {output} {input}"

rule counts:
    input:
        tra = "output/transcripts.idx",
        fq1 = lambda wildcards: table['fq1'][wildcards.sample],
        fq2 = lambda wildcards: table['fq2'][wildcards.sample]
    output:
        directory("output/abundances/{sample}")
    shell:
        "kallisto quant -i {input.tra} -o {output} -b 100 {input.fq1} {input.fq2}"

# rule normalize:
#    input:
#           expand("output/abundances/{sample}/abundace.h5", sample=SAMPLES)
#    output:
# mal gucken
#    script:
#        "scripts/sleuth.R"

#rule boxenplot:
