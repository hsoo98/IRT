import json
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_table('../DEG_table_ttest.tsv')
with open ('../Marker_Dict.json') as f:
	Marker_Dict = json.load(f)

SIGNIFICANT = True # False

if SIGNIFICANT:
	#data = data.loc[(data.Pvalue<0.05)&(data.Log2FC.abs()>=1)]
	data = data.loc[data.Pvalue<0.05]

for marker in Marker_Dict:
	plt.figure(figsize=(8,8))
	sns.scatterplot(x='WT_avg', y='CEB_KO_avg', data=data, c='gray')
	# Gray dots only for the significant genes?

	Tfr_rel = data.loc[data.GeneName.isin(Marker_Dict[marker])]

	Tfr_pos = Tfr_rel.loc[Tfr_rel.WT_avg<Tfr_rel.CEB_KO_avg]
	Tfr_neg = Tfr_rel.loc[Tfr_rel.WT_avg>Tfr_rel.CEB_KO_avg]

	sns.scatterplot(x='WT_avg', y='CEB_KO_avg', data=Tfr_neg, c='blue')
	sns.scatterplot(x='WT_avg', y='CEB_KO_avg', data=Tfr_pos, c='red')

	plt.xscale('log')
	plt.yscale('log')

	for i in Tfr_pos.index:
		plt.text(x=Tfr_pos.WT_avg[i]+1, y=Tfr_pos.CEB_KO_avg[i]+1, s=Tfr_pos.GeneName[i],
		         fontdict=dict(color='red', size=10))
	for i in Tfr_neg.index:
		plt.text(x=Tfr_neg.WT_avg[i]+1, y=Tfr_neg.CEB_KO_avg[i]+1, s=Tfr_neg.GeneName[i],
	    	     fontdict=dict(color='blue', size=10))
	if SIGNIFICANT:
		plt.savefig('../output/Fig1_Scatter_ttest_%s_sig.pdf'%marker)
	else:
		plt.savefig('../output/Fig1_Scatter_ttest_%s.pdf'%marker)
	plt.show()
