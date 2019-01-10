import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sleuth_table = pd.read_table("temp/sleuth_table.tsv", delim_whitespace=True)

sleuth_table = sleuth_table.sort_values('pval', ascending=True)
sleuth_table = sleuth_table.head(20)



p_value_table = sleuth_table[['target_id', 'pval']]
print(p_value_table)
