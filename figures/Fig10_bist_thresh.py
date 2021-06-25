import numpy as np
import matplotlib.pyplot as plt
from sbcf import *
import pickle

#membership
mmax = 4
gm = np.zeros(mmax+1)
gm[mmax] += 1

#infection
beta = lambda n,i,trate,nu: trate*i**nu

# nmax_list = [40]
nmax_list = [n for n in range(3,41)] #overflow problem beyon 41
gamma_list = np.linspace(2.01,5.,51)
#add some more value for fun
plot_gamma = [2.1,3.,3.8,4.,4.2]
gamma_list = np.concatenate((gamma_list,plot_gamma))
gamma_list.sort()
results = {nmax:dict() for nmax in nmax_list}


for nmax in nmax_list:
    print("----------")
    print(f"nmax: {nmax}")
    print("----------")
    #heterogeneous groups
    n = np.arange(nmax+1)
    nmin = 2

    for i,gamma in enumerate(gamma_list):
        print(f"gamma: {gamma}")
        print("--------------")
        pn = np.zeros(nmax+1)
        pn[nmin:] = n[nmin:]**(-gamma)
        pn[0:nmin] = 0.
        pn /= np.sum(pn)

        nu_c = bistability_threshold_safe(beta, gm, pn, max_params=(1,4.5))

        print("bistability threshold: ",nu_c)
        lambda_c = invasion_threshold_safe(beta, gm, pn, fixed_args=(nu_c,))
        print("invasion threshold: ",lambda_c)

        results[nmax][gamma] = nu_c

    # # plt.plot(gamma_list,[results[nmax][gamma] for gamma in gamma_list],
            # # label=r"$m_\mathrm{max}=$"+f"${nmax}$")
# for gamma in plot_gamma:
    # plt.semilogx(nmax_list,[results[nmax][gamma] for nmax in nmax_list],
            # label=fr"$\gamma={gamma}$")
# plt.legend()
# plt.ylabel(r'Bistability threshold $\nu_\mathrm{c}$')
# # plt.xlabel(r'Membership distribution exponent $\gamma_m$')
# plt.xlabel(r'Maximal membership $m_\mathrm{max}$')
# plt.show()

with open('./dat/bistability_heterogeneous_groups.pk', 'wb') as outfile:
    pickle.dump(results,outfile)
