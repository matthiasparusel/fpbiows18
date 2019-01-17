import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sleuth_table = pd.read_table(snakemake.input[0], delim_whitespace=True)
p_value_table = sleuth_table['pval']
ranges = [round(0.01*x, 2) for x in range(0,105,5)] # 5% intervalls

sns.set(style='darkgrid')
ax = sns.distplot(p_value_table, bins=ranges, axlabel='P-values', hist=True, rug=False, kde=False)
ax.yaxis.grid(True)
ax.xaxis.grid(True)
plt.ylabel('counts')
plt.savefig(snakemake.output[0])
