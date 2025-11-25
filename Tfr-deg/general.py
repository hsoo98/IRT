import pandas as pd
import numpy as np
from tqdm import tqdm
import json

import math
from scipy.stats import mannwhitneyu as mwu
from scipy.stats import ttest_ind

Sample_Dict = {
	'WT': ['D21_WT_Ttr_129K','D21_WT_Ttr_196K','D21_WT_Ttr_275K'],
	'CEB_KO':['D14_CEB_KD_192K','D21_CEB_KD_123K']
}

Marker_Dict = {
	'Tfr_functional':[
		'Foxp3','Bcl6','Cxcr5','Ctla4','Icos','Tnfrsf18','Pdcd1','Prdm1','Btla','Cd83',
		'Bach2','Tnfrsf9','Tnfrsf1b','Tgfb1'
	],
	'Tfh_related':[
		'Tox2','Tox','Il6st','Egr2','Sh2d1a','Foxo1','Pou2af1','Il21r','Peli2','Il17ra',
		'Klf2','Tgfbr2','Bach2','Ikzf3','Id3'
	],
	'Treg_related':[
		'Gzmb','Ikzf2','Prdm1','Il2rb','Il2ra','Tnfrsf9'
	],
	'Naive_Treg':[
		'Klf2','S1pr1','S1pr4','Foxo1','Sell'
	],
	'Th1_related':[
		'Tbx21','Ifng','Cxcr3','Il23r','Il12rb1','Il12rb2','Cxcr6'
	],
	'TCR_signaling':[
		'Cd3g','Cd3d','Cd3e','Lck'
	],
	'Activation':[
		'Mki67'
	],
	'Th17_related':[
		'Cxcr4','Rorc'
	],
	'Autoimmune':[
		'Ly6a'
	],
	'Cytokines':[
		'Il16','Il6','Il1b','Il21','Il2','Il12a','Il15','Il18','Il7','Il25',
		'Il23r','Il17d','Il18bp','Il22b'
	]
}

with open ('../Marker_Dict.json','w') as f:
	json.dump(Marker_Dict, f)

def GetDatasets():
	Datasets = {}
	Genes = set([])
	for group in ['WT', 'CEB_KO']:
		for sample in Sample_Dict[group]:
			df = pd.read_table('../StringTie/%s_gene_abundances.tsv'%sample)
			df = df.loc[~df['Gene Name'].isin(['-','.'])]
			df = df.loc[~((df['Gene Name'].str[:2]=='Gm')&
			              (df['Gene Name'].str[3].isin(list(map(lambda x: str(x),range(10))))))]
			df = df.loc[df.Coverage>1]
			Datasets[sample] = df
			Genes = Genes.union(set(list(df['Gene Name'])))
	print("%d genes obtained"%len(Genes))
	return Datasets, Genes

def DEGcalc(method='ttest'):
	# DESeq2 is only for the count data --> check the difference
	# MannWhitneyU test & ttest_rel both show non-significance
	Datasets, Genes = GetDatasets()
	output = {'GeneName':[], 'WT_avg':[], 'CEB_KO_avg':[], 'Log2FC':[], 'Pvalue':[]}
	for gene in tqdm(Genes):
		wt = []
		for sample in Sample_Dict['WT']:
			df = Datasets[sample]
			if gene in list(df['Gene Name']):
				wt.append(df.loc[df['Gene Name']==gene]['TPM'].tolist()[0])
		cebko = []
		for sample in Sample_Dict['CEB_KO']:
			df = Datasets[sample]
			if gene in list(df['Gene Name']):
				cebko.append(df.loc[df['Gene Name']==gene]['TPM'].tolist()[0])
		if len(wt)==0 or len(cebko)==0: continue
		if method=='ttest':
			_, pval = ttest_ind(wt, cebko)
		else:
			_, pval = mwu(wt, cebko)
		l2fc = math.log2(np.mean(cebko)/np.mean(wt))
		output['GeneName'].append(gene)
		output['WT_avg'].append(np.mean(wt))
		output['CEB_KO_avg'].append(np.mean(cebko))
		output['Log2FC'].append(l2fc)
		output['Pvalue'].append(pval)
	output = pd.DataFrame(output)
	print('%d genes were remained'%len(output))
	return output, method

if __name__ == "__main__":
	output, method = DEGcalc()
	output.to_csv('../DEG_table_%s.tsv'%method,sep='\t',index=False)
