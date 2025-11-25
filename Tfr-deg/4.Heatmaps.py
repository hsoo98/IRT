import json
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

from general import GetDatasets, Marker_Dict, Sample_Dict

Datasets, Genes = GetDatasets()

for marker in Marker_Dict:
	output = {'Gene':Marker_Dict[marker]}
	for group in Sample_Dict:
		for i, sample in enumerate(Sample_Dict[group]):
			prefix = '%s%d'%(group,i)
			exp_list = []
			df = Datasets[sample]
			for gene in Marker_Dict[marker]:
				if gene not in list(df['Gene Name']):
					exp_list.append(np.nan)
				else:
					exp_list.append(df.loc[df['Gene Name']==gene]['TPM'].tolist()[0])
			output[prefix] = exp_list # Normalization??
	output = pd.DataFrame(output).dropna()
	output = output.set_index('Gene', drop=True)

	sns.heatmap(output, cmap='bwr')
	plt.savefig('../Fig4_Heatmap_%s.pdf'%marker)
	plt.show()
