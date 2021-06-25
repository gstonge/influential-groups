import pickle
import numpy as np
import matplotlib.pyplot as plt

# datafile = './socio_data/form_thr3_LyonSchool.pk'
datafile = './socio_data/form_coauth-dblp_sub.pk'

with open(datafile,'rb') as filename:
    data = pickle.load(filename)

m_list = data['m_list']
n_list = data['n_list']
nmax = data['nmax']
mmax = data['mmax']
gm = data['gm']
pn = data['pn']

print(len(m_list))
print(len(n_list))
print(np.mean(m_list))
print(np.mean(n_list))
print(mmax)
print(nmax)


plt.scatter(np.arange(mmax+1)[1:],gm[1:])
print(gm[-1])
plt.xscale('log')
plt.yscale('log')
plt.ylim([10**-5,1])
plt.show()

plt.bar(np.arange(nmax+1),pn)
plt.show()
