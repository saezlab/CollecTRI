# The CollecTRI GRN is implemented in pypath/Omnipath (using the CollecTRI_signed.csv).
# Since pypath includes some filtering and quality control steps
# e.g. checking conversion of gene names to uniprot ids
# we load CollecTRI through pypath and continue our analysis with the final
# CollecTRI GRN

from pypath.core import network
from pypath.resources import network as netres
from pypath.omnipath import export
import pandas as pd
import numpy as np

n = network.Network(netres.collectri, allow_loops = True)
n.load(netres.tf_mirna['collectri'])
e = export.Export(n)
e.make_df(unique_pairs=False)

e.df['weight'] = np.where(e.df['consensus_stimulation'] == 1, 1, np.where(e.df['consensus_inhibition'] == 1, -1, np.nan))
e.df = e.df.loc[:, ['source_genesymbol', 'target_genesymbol', 'weight', 'sources', 'references']].rename(columns={'source_genesymbol': 'source', 'target_genesymbol': 'target', 'sources': 'resources'})

pd.DataFrame.to_csv(e.df, "output/CollecTRI/CollecTRI_GRN.csv", index = False)
