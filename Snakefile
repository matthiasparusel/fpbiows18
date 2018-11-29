rule boxenplot:
    input:
        "data/reads/{sample}.fq"
    output:
        "read_counts/{sample}.svg"
    script:
        "scripts/boxenblots.py"
