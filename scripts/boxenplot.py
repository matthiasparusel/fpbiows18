# import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# import sys

# Get the sample read
samples = pd.read_table(snakemake.input[0], delim_whitespace=True)

# First, get Kallisto results and estimated counts
kal_results = list()
for sample in samples['sample']:
    df = pd.read_table(f'./results/{sample}/kallisto/abundance.tsv',
                       delim_whitespace=True)
    df['sample'] = sample
    kal_results.append(df)
kal_results = pd.concat(kal_results)

# Second, get the normalized estimated counts
sle_results = pd.read_table(snakemake.input[1], delim_whitespace=True)

# Finally plot the thing
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 4), sharey=True)
sns.boxenplot(x='sample', y='est_counts', data=kal_results, ax=axes[0])
axes[0].set_title('counts of all samples')
sns.boxenplot(x='sample', y='est_counts', data=sle_results, ax=axes[1])
axes[1].set_title('normalized counts of all samples')
for ax in axes:
    ax.yaxis.grid(True)
plt.savefig(snakemake.output[0])
