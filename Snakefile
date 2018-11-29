rule boxenplot:
    input:
        "data/table.tsv",
        "data/reads/{sample}.1.fq",
        "data/reads/{sample}.2.fq"
    output:
        "read_counts/{sample}.svg"
    script:
        "scripts/boxenblots.py"
