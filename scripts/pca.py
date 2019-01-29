from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sample_table = pd.read_table(snakemake.input[0], delim_whitespace=True)
sleuth_matrix = pd.read_table(snakemake.input[1], delim_whitespace=True)

sleuth_matrix = sleuth_matrix.transpose()
pca = PCA(n_components=2)
X = pca.fit(sleuth_matrix)
X = X.transform(sleuth_matrix)

pca_df = pd.concat([pd.DataFrame(X), sample_table[['sample', 'condition']]], axis=1)

sns.set(style='darkgrid')
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), sharey=True)
sns.scatterplot(x=0, y=1, data=pca_df, legend='full', hue='condition', ax=axes[0])
sns.scatterplot(x=0, y=1, data=pca_df, legend='full', hue='sample', ax=axes[1])
for ax in axes:
    ax.set_xlabel('')
    ax.set_ylabel('')

plt.savefig(snakemake.output[0], bbox_inches="tight")
