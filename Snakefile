import pandas as pd

table = pd.read_table("data/table.tsv").set_index("sample", drop=False)
SAMPLES = table.index.values.tolist()

rule all:
    input:
        expand("output/abundances/{sample}", sample=SAMPLES)


rule index:
    input:
        "data/ref/transcriptome.chr21.fa"
    output:
        "output/transcripts.idx"
    shell:
        "kallisto index -i {output} {input}"


rule counts:
    input:
        tra = "output/transcripts.idx", # can stay
        read1 = lambda wildcards: table['fq1'][wildcards.sample],
        read2 = lambda wildcards: table['fq2'][wildcards.sample]
    output:
        directory("output/abundances/{sample}")
    shell:
        "kallisto quant -i {input.tra} -o {output} -b 100 {input.read1} {input.read2}"


# rule normalize:
    input:
        expand("output/abundances/{sample}/abundace.h5", sample=SAMPLES)
#    output:
# mal gucken
#    script:
        "scripts/sleuth.R"

#rule boxenplot:
