import numpy as np
import matplotlib.pyplot as plt
from sbcf import *
import pickle

#structure
mmax = 4
gm = np.zeros(mmax+1)
gm[mmax] += 1

nmax = 20
nmin = 2
n= np.arange(nmax+1)
pn = np.zeros(nmax+1)

#infection
beta = lambda n,i,trate,nu: trate*i**nu

nu_list = [1., 1.5, 2.]
gamma_list = np.linspace(2.01,5.5,101)
#add some more value for fun
results = {nu:dict() for nu in nu_list}

for nu in nu_list:
    print(f'nu:{nu}')
    for gamma in gamma_list:
        pn[nmin:] = n[nmin:]**(-gamma)
        pn /= np.sum(pn)
        results[nu][gamma] = invasion_threshold(
            beta, gm, pn,fixed_args=(nu,), initial_param=0.2)

    plt.plot(gamma_list,[results[nu][gamma] for gamma in gamma_list],
            label=r"$\nu=$"+f"${nu}$")
plt.legend()
plt.show()

with open('./dat/Fig6_exp.pk', 'wb') as outfile:
    pickle.dump(results,outfile)
