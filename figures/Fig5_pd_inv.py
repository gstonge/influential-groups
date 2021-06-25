import numpy as np
import matplotlib.pyplot as plt
from gcm import *
import pickle
from scipy.special import loggamma
from scipy.optimize import fsolve



#membership
mmax = 3
gm = np.zeros(mmax+1)
gm[mmax] += 1


#group distribution
nmax = 4
pn = np.zeros(nmax+1)
pn[nmax] += 1

#infection
beta = lambda n,i,trate,nu: trate*i**nu

#get the tricritical point
nu_c = bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                              max_params=(1,7))
lambda_c = invasion_threshold_safe(beta, gm, pn, fixed_args=(nu_c,),
                            min_param=10**(-14), max_param=1)
tricrit = (lambda_c,nu_c)
print(f"Tricritical point : {tricrit}")

#get the invasion threshold in a line, passing through tricritical points
nu_list = np.linspace(0.1,4.,100)

lambda_list = []
for nu in nu_list:
    lambda_c = invasion_threshold_safe(beta, gm, pn, fixed_args=(nu,),
                            min_param=10**(-14), max_param=1)
    lambda_list.append(lambda_c)

results = dict()
results['lambda_list'] = lambda_list
results['nu_list'] = nu_list
with open('./dat/Fig3_pd_inv.pk', 'wb') as filename:
    pickle.dump(results,filename)

plt.plot(lambda_list,nu_list)
plt.scatter([tricrit[0]],[tricrit[1]], marker='*')
plt.show()
