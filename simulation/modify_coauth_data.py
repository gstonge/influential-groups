"""
Construct hypergraphs with the same group size sequence as the coauthorship
data, but with a homogeneous membership distribution
"""

import pickle
import numpy as np
import horgg
import matplotlib.pyplot as plt
from scipy.special import loggamma


with open('./socio_data/form_coauth-dblp_sub.pk', 'rb') as filename:
    data = pickle.load(filename)


m_list = data['m_list']
n_list = data['n_list']
pn = data['pn']
nmax = max(n_list)
mmean = np.mean(m_list)
mmax = 4
gm = np.zeros(mmax+1)
gm[4] = mmean-3
gm[3] = 1-gm[4]

#new values
seed = 42
np.random.seed(seed)

m_list = horgg.sequence_2(n_list,gm)
mmax = max(m_list)
gm = np.zeros(mmax+1)
for m in m_list:
    gm[m] += 1
gm /= np.sum(gm)

#save new modified dataset
results = dict()
results['m_list'] = m_list
results['n_list'] = n_list
results['gm'] = gm
results['pn'] = pn
results['nmax'] = nmax
results['mmax'] = mmax
results['edge_list'] = [] #dummy

with open('./socio_data/form_modified_coauth-dblp_sub.pk', 'wb') as filename:
    pickle.dump(results,filename)

