import pandas as pd
import numpy as np
import os
import sys

from scipy.stats import mannwhitneyu as mwu
from scipy.stats import ttest_ind
import math

import matplotlib.pyplot as plt
import seaborn as sns

gse = pd.read_table('../GSE124883/GSE124883_Tissue_Matrix.txt')

try: cutoff = float(sys.argv[1])
except: cutoff = 0.001
print('P-value cutoff: %f'%cutoff)

# LN data --> for the volcano plot
# SP data --> for obtaining the gene sets
# Tcon vs Tfr
# Tfh vs Tfr
# Treg vs Tfr

def DEGcalc (df, tissue, target1, target2, method='ttest'):
	output = {'GeneID':[], '%s_exp'%target1:[], '%s_exp'%target2:[], 'Log2FC':[], 'Pvalue':[]}
	for i, row in df.iterrows():
		exp1 = row[list('%s %s %d'%(tissue, target1, x) for x in range(1,4))].tolist()
		exp2 = row[list('%s %s %d'%(tissue, target2, x) for x in range(1,4))].tolist()
		exp1_avg = np.mean(exp1)
		exp2_avg = np.mean(exp2)
		if exp1_avg == 0 or exp2_avg ==0: continue
		if method == 'ttest':
			_, pval = ttest_ind(exp1, exp2)
		else:
			print(method)
			_, pval = mwu(exp1, exp2)
		l2fc = math.log2(exp2_avg/exp1_avg)
		output['GeneID'].append(row['id'])
		output['%s_exp'%target1].append(exp1_avg)
		output['%s_exp'%target2].append(exp2_avg)
		output['Log2FC'].append(l2fc)
		output['Pvalue'].append(pval)
	output = pd.DataFrame(output)
	output.to_csv('../GSE124883/GSE124883_DEG_%s_%s_vs_%s.tsv'%(tissue, target1, target2),
	              sep='\t',index=False)
	return output

def Volcanoplot (df, genelist, t1, t2):
	df['-logP'] = df['Pvalue'].apply(lambda x: -math.log10(x))
	plt.figure(figsize=(8,8))
	sns.scatterplot(x='Log2FC', y='-logP', data=df, c='lightgray')
	
	try: target = df.loc[df.GeneID.isin(genelist)]
	except: target = df.loc[df.GeneName.isin(genelist)]
	sns.scatterplot(x='Log2FC', y='-logP', data=target, c='red')

	sig_target = target.loc[target.Pvalue<cutoff]

	#for i in sig_target.index:
	#	plt.text(x=sig_target['Log2FC']+1, y=sig_target['-logP']+1, s=sig_target.GeneID[i],
	#	         fontdict=dict(color='red', size=10))
	plt.savefig('../output/Fig2_Volcano_%s_%s_%.3f.pdf'%(t1, t2, cutoff))
	#plt.show()
	return sig_target

for match in [('Tcon', 'Tfr'), ('Tfh', 'Tfr'), ('Treg', 'Tfr')]:
	target1 = match[0]
	target2 = match[1]
	if 'GSE124883_DEG_Spleen_%s_vs_%s.tsv'%(target1, target2) in os.listdir('../GSE124883'):
		df1 = pd.read_table('../GSE124883/GSE124883_DEG_Spleen_%s_vs_%s.tsv'%(target1, target2))
	else:
		df1 = DEGcalc(gse, 'Spleen', target1, target2)
	genelist = df1.loc[(df1.Log2FC>0)&(df1.Pvalue<cutoff)].GeneID.tolist()
	print('%d DEGs between %s vs %s'%(len(genelist), target1, target2))
	
	#if 'GSE124883_DEG_LN_%s_vs_%s.tsv'%(target1, target2) in os.listdir('../GSE124883'):
	#	df2 = pd.read_table('../GSE124883/GSE124883_DEG_LN_%s_vs_%s.tsv'%(target1, target2))
	#else:
	#	df2 = DEGcalc(gse, 'LN', target1, target2)
	df2 = pd.read_table('../DEG_table_ttest.tsv')
	output = Volcanoplot (df2, genelist, target1, target2)
	output.to_csv('../output/Fig2_Targets_%s_%s_%.3f.tsv'%(target1, target2, cutoff),
	              sep='\t', index=False)

