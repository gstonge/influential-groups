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

state_meta = get_state_meta(mmax, nmax, gm, pn)

#infection
beta = lambda n,i,trate,nu: trate*i**nu

#get the tricritical point
nu_c = bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                              max_params=(1,7))
lambda_c = invasion_threshold_safe(beta, gm, pn, fixed_args=(nu_c,),
                            min_param=10**(-14), max_param=1)
tricrit = (lambda_c,nu_c)
print(f"Tricritical point : {tricrit}")

#get the persistence threshold in a line, passing through tricritical points
nu_list = np.linspace(0.1,4.,100)

lambda_list = []
for nu in nu_list:
    print(f"---------")
    print(f"nu: {nu}")
    lambda_c = invasion_threshold_safe(beta, gm, pn, fixed_args=(nu,),
                            min_param=10**(-14), max_param=1)
    #start from invasion threshold
    param_init = 1.1*lambda_c
    param_var = -0.001*lambda_c
    param_list,stationary_state_list,infected_fraction_list = stable_branch(
        beta, state_meta, param_init, param_var, fixed_args=(nu,),
        h=0.01, rtol=10**(-10), Jtol=10**(-4),
        min_density=10**(-6),max_iter=1000, verbose=False)
    print(f"I: {infected_fraction_list[-1]}")
    print(f"lambda: {param_list[-1]}")
    print(f"---------")
    lambda_list.append(param_list[-1])

results = dict()
results['lambda_list'] = lambda_list
results['nu_list'] = nu_list
with open('./dat/Fig3_pd_pers.pk', 'wb') as filename:
    pickle.dump(results,filename)

plt.plot(lambda_list,nu_list, '--')
plt.scatter([tricrit[0]],[tricrit[1]], marker='*')
plt.show()
