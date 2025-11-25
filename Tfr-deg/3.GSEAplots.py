import pandas as pd
import numpy as np
import sys

import matplotlib.pyplot as plt
from gseapy.plot import gseaplot
import gseapy as gp

try: cutoff = float(sys.argv[2])
except: cutoff = 0.01
print('Pvalue < %.3f'%cutoff)

# Prefix 1
database = 'DisGeNET'
target_term = 'Lupus Erythematosus, Systemic' 
# 'Lupus Erythematosus', 'Lupus Erythematosus, Discoid', 'Lupus Vulgaris'

# Prefix 2
database = 'KEGG_2021_Human'
target_term = 'Systemic lupus erythematosus'

def GSEAPYplot (df, database, target_term):
	df = df.loc[df.Pvalue<cutoff]
	df['Rank'] = -np.log10(df.Pvalue)*df.Log2FC
	df = df.sort_values('Rank', ascending=False).reset_index(drop=True)
	df.GeneName = df.GeneName.str.upper()  # May need better ortholog mapping 
	ranking = df[['GeneName','Rank']]

	pre_res = gp.prerank(rnk=ranking, gene_sets=database, seed=6)
	out = []
	for term in list(pre_res.results):
		out.append([term,
		            pre_res.results[term]['fdr'],
					pre_res.results[term]['es'],
					pre_res.results[term]['nes']])
	out_df = pd.DataFrame(out, columns=['Term','fdr','es','nes']).sort_values('fdr')
	out_df = out_df.reset_index(drop=True)
	gseaplot(term = target_term, **pre_res.results[target_term])
	plt.savefig('../Fig3_GSEA_(%s)_(%s)_%.3f.pdf'%(database, target_term, cutoff))
	plt.show()
	return out_df

def GSEAPYplot_custom (df, prefix, genelist):
	df = df.loc[df.Pvalue<cutoff]
	df['Rank'] = -np.log10(df.Pvalue)*df.Log2FC
	df = df.sort_values('Rank', ascending=False).reset_index(drop=True)
	df.GeneName = df.GeneName.str.upper()
	ranking = df[['GeneName','Rank']]
	
	user_set = {prefix: genelist}
	pre_res = gp.prerank(rnk=ranking, gene_sets=user_set, seed=6)
	gseaplot(term=prefix, **pre_res.results['custom']
	plt.savefig('../Fig3_GSEA_(custom)_(%s)_%.3f.pdf'%(prefix, cutoff))
	plt.show()

if __name__ == "__main__":
	df = pd.read_table('../DEG_table_ttest.tsv')
	GSEAPYplot(df, database, target_term)
