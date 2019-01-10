import seaborn as sns
import panda as pd
import matplotlib.pyplot as plt

sample = pd.read_table(snakemake.input[0])

sns.boxenplot(data = sample)
plt.figure()
plt.savefig(snakemake.output.counts)
