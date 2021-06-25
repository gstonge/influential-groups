import numpy as np
import matplotlib.pyplot as plt
from sbcf import *
import pickle

#structure
mmax = 4
# mmax = 10
gm = np.zeros(mmax+1)
gm[mmax] += 1

nmax = 20
nmin = 2

#infection
beta = lambda n,i,trate,nu: trate*i**nu

nu_list = [1., 1.5, 2.]
gamma = 3.
nmax_list = [i for i in range(10,101)]
#add some more value for fun
results = {nu:dict() for nu in nu_list}

for nu in nu_list:
    print(f'nu:{nu}')
    for nmax in nmax_list:
        n= np.arange(nmax+1)
        pn = np.zeros(nmax+1)
        pn[nmin:] = n[nmin:]**(-gamma)
        pn /= np.sum(pn)
        results[nu][nmax] = invasion_threshold_safe(
            beta, gm, pn,fixed_args=(nu,), max_param=0.1)

    plt.loglog(nmax_list,[results[nu][nmax] for nmax in nmax_list],
            label=r"$\nu=$"+f"${nu}$")
plt.loglog(nmax_list, 2*np.array(nmax_list,dtype=np.float64)**(-1.5), '--')
plt.legend()
plt.show()

# with open('./dat/Fig6_size.pk', 'wb') as outfile:
    # pickle.dump(results,outfile)
