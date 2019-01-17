import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sleuth_table = pd.read_table("temp/sleuth_table.tsv", delim_whitespace=True)
sleuth_table = sleuth_table.sort_values('pval', ascending=True)
sleuth_table = sleuth_table.head(20)

kallisto_table = pd.read_table("temp/counts_normalized.tsv", delim_whitespace=True)

p_value_table = sleuth_table[['target_id', 'pval']]
counts_table = kallisto_table[['target_id', 'est_counts', 'sample']]



merged_table = counts_table.merge(p_value_table, left_on='target_id', right_on='target_id', how='inner')


plt.figure(figsize=[12,8])
ax = sns.stripplot(x='est_counts', y='pval', hue='target_id', data=merged_table)

plt.xticks(rotation=45)
plt.grid()
plt.tight_layout()
plt.savefig('results/stripplot.pdf')

#print(p_value_table)

#print(merged_table)
