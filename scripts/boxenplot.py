import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

#Argumente einlesen (Format: python boxenplot.py pfad1 pfad2 ... pfadn)
paths = sys.argv


kallisto_output_a =  pd.read_table('output/a_abundance/abundance.tsv', delim_whitespace=True)
kallisto_output_b =  pd.read_table('output/b_abundance/abundance.tsv', delim_whitespace=True)

counts_a = kallisto_output_a.loc[: ,"est_counts"]
counts_b = kallisto_output_b.loc[: ,"est_counts"]

ax = sns.boxenplot(data=[counts_a, counts_b])

plt.show()
