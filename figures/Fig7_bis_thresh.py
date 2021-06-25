import numpy as np
import matplotlib.pyplot as plt
from sbcf import *
import pickle

#group
nmax = 4
pn = np.zeros(nmax+1)
pn[nmax] += 1

#infection
beta = lambda n,i,trate,nu: trate*i**nu

mmax_list = [100,300,1000,3000,10000]
gamma_list = np.linspace(2.01,5.5,101)
#add some more value for plot
plot_gamma = [2.1,3.,3.8,4.,4.2]
gamma_list = np.concatenate((gamma_list,plot_gamma))
gamma_list.sort()
results = {mmax:dict() for mmax in mmax_list}

for mmax in mmax_list:
    #heterogeneous membership
    mmin = 2
    m = np.arange(mmax+1)

    for gamma in gamma_list:
        # print(f"gamma: {gamma}")
        # print("--------------")
        gm = np.zeros(mmax+1)
        gm[mmin:] = m[mmin:]**(-gamma)
        gm[0:mmin] = 0.
        gm /= np.sum(gm)
        nu_c = bistability_threshold(beta, gm, pn, initial_params=(0.1,4.8))
        # print("bistability threshold: ",nu_c)

        results[mmax][gamma] = nu_c


with open('./dat/Fig5_bis_thresh.pk', 'wb') as outfile:
    pickle.dump(results,outfile)
