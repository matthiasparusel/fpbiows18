import pandas as pd

table = pd.read_table("data/table.tsv").set_index("sample", drop=False)
# fq_a = table['fq1']['a'] + ' ' + table['fq2']['a']
# fq_b = table['fq1']['b'] + ' ' + table['fq2']['b']

rule index:
    input:
        "data/ref/transcriptome.chr21.fa"
    output:
        "output/transcripts.idx"
    shell:
        "kallisto index -i {output} {input}"

rule counts:
    input:
        tra="output/transcripts.idx",
        a1="data/reads/a.chr21.1.fq",
        a2="data/reads/a.chr21.2.fq"
    output:
        directory("output/a_abundance")
    shell:
        "kallisto quant -i {input.tra} -o {output} -b 100 {input.a1} {input.a2}"
