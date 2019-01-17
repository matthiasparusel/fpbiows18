import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sleuth_table = pd.read_table(snakemake.input[0], delim_whitespace=True)
sleuth_table = sleuth_table.sort_values('pval', ascending=True)
sleuth_table = sleuth_table.head(20)

kallisto_table = pd.read_table(snakemake.input[1], delim_whitespace=True)

p_value_table = sleuth_table[['target_id', 'pval']]
counts_table = kallisto_table[['target_id', 'est_counts', 'sample']]



merged_table = counts_table.merge(p_value_table, left_on='target_id', right_on='target_id', how='inner')


plt.figure(figsize=[12,8])
ax = sns.stripplot(x='est_counts', y='pval', hue='target_id', data=merged_table)

plt.xticks(rotation=45)
plt.grid()
plt.tight_layout()
plt.savefig(snakemake.output[0])

#print(p_value_table)

#print(merged_table)
