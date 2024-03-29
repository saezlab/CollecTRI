# The CollecTRI GRN is implemented in pypath/Omnipath (using the CollecTRI_signed.csv).
# Since pypath includes some filtering and quality control steps
# e.g. checking conversion of gene names to uniprot ids
# we load CollecTRI through pypath and continue our analysis with the final
# CollecTRI GRN

from pypath.core import network
from pypath.omnipath import export
from pypath.share import common
from pypath.resources import network as netres
from pypath.internals import intera
import pandas as pd
import numpy as np

n = network.Network(netres.collectri, allow_loops = True)
n.load(netres.tf_mirna['collectri'])
e = export.Export(n)
e.make_df(unique_pairs=False, extra_edge_attrs = {'sign_decision': lambda e, d: ';'.join({sd for ev in e.direction[d] for sd in common.to_set(ev.attrs.get('sign_decision'))})})
indx_int = [not value for value in [isinstance(s, intera.Complex) for s in e.df.source]]
indx_comp = [isinstance(s, intera.Complex) for s in e.df.source]

e.df['weight'] = np.where(e.df['is_stimulation'] == 1, 1, np.where(e.df['is_inhibition'] == 1, -1, np.nan))
e.df = e.df.loc[:, ['source_genesymbol', 'target_genesymbol', 'weight', 'sources', 'references', 'sign_decision']].rename(columns={'source_genesymbol': 'source', 'target_genesymbol': 'target', 'sources': 'resources'})

#split data frame into interactions from complexes to merge them again
intdf = e.df[indx_int].copy()
compdf = e.df[indx_comp].copy()
compdf['source'] = compdf['source'].astype(str)

#return to complex name
#AP1
compdf.loc[compdf['source'].str.startswith('JUN') | compdf['source'].str.startswith('FOS'), 'source'] = 'AP1'

#NFKB
compdf.loc[compdf['source'].str.startswith('REL') | compdf['source'].str.startswith('NFKB'), 'source'] = 'NFKB'

full_df = pd.concat([intdf, compdf], axis=0)
full_df['sign_decision'] = full_df['sign_decision'].astype(str)
full_df.drop_duplicates(inplace=True)
full_df['sign_decision'] = full_df['sign_decision'].str.replace('PMID;default activation', 'PMID', regex=False) #for NFKB1:DEFB4A due to gene name conversions
pd.DataFrame.to_csv(full_df, "output/CollecTRI/CollecTRI_GRN.csv", index = False)
