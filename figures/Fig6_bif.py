import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import loggamma
from gcm import *

compute = True

#infection
nu = 2.80
beta = lambda n,i,trate,nu: trate*i**nu

#Poisson membership
mmin = 1
nmax = 3
mmax = 100
m = np.arange(mmax+1)
pn = np.zeros(nmax+1)
pn[nmax] += 1
param = 5
gm = np.exp(m*np.log(param) - loggamma(m+1))
gm[0:mmin] = 0.
gm /= np.sum(gm)
state_meta_1 = get_state_meta(mmax, nmax, gm, pn)
assert np.isclose(np.sum(m*(m-1)*gm)/np.sum(m*gm),5)
print(np.sum(gm*m))
print(np.sum(gm*m**2))
print(np.sum(gm*m**3))
print("invasion threshold: ",
      invasion_threshold(beta, gm, pn,fixed_args=(nu,), initial_param=0.2))
print("bistability threshold: ",
      bistability_threshold(beta, gm, pn, initial_params=(0.1,1.1)))

#heterogeneous membership
mmin = 4
nmax = 3
mmax = 100
m = np.arange(mmax+1)
pn = np.zeros(nmax+1)
pn[nmax] += 1

a = -0.8
b = 3.72
c = 40
gm = np.zeros(mmax+1)
gm[mmin:] = np.exp(-m[mmin:]/c)/(a+m[mmin:])**b
gm[0:mmin] = 0.
gm /= np.sum(gm)
state_meta_2 = get_state_meta(mmax, nmax, gm, pn)
# assert np.isclose(np.sum(m*(m-1)*gm)/np.sum(m*gm),5)
print(np.sum(gm*m))
print(np.sum(gm*m**2))
print(np.sum(gm*m**3))
print("invasion threshold: ",
      invasion_threshold(beta, gm, pn,fixed_args=(nu,), initial_param=0.2))
print("bistability threshold: ",
    bistability_threshold(beta, gm, pn, initial_params=(0.1,1.1)))

state_meta_dict = {"struct1":state_meta_1,
                   "struct2":state_meta_2}

if compute:
    #prepare result dict
    results = {key:dict() for key in state_meta_dict}
    for key,state_meta in state_meta_dict.items():
        results[key]['stable'] = dict()
        results[key]['gm'] = list(state_meta[3])
        if key == "struct1":
            results[key]['unstable'] = dict()
    #compute stable and unstable solutions
    for key,state_meta in state_meta_dict.items():
        print(key)
        gm = state_meta[3]
        pn = state_meta[4]
        param_c = invasion_threshold(beta,gm,pn,fixed_args=(nu,))
        param_init = 1.1*param_c

        #get the stable branch
        param_var = -0.005 #backward direction
        Jtol = 0.0001
        stable_param,stable_fixed_point,stable_infected_fraction = \
                stable_branch(beta,state_meta,param_init,param_var,
                              rtol=10**(-12), max_iter=20000,
                              fixed_args=(nu,),Jtol=Jtol,verbose=False)
        stable_param = [param/param_c for param in stable_param]

        results[key]['stable']['param'] = stable_param
        results[key]['stable']['infected_fraction'] = stable_infected_fraction
        if key == "struct1":
            #get unstable branch
            param_var = abs(stable_param[-2]-stable_param[-1])*param_c
            fni = stable_fixed_point[-1][1]
            param_init = stable_param[-1]*param_c
            unstable_param,unstable_fixed_point,unstable_infected_fraction = \
                    unstable_branch(fni,beta,state_meta,param_init,
                                    param_var,fixed_args=(nu,),init_iter=100,
                                    h=10**(-2),
                                    rtol=10**(-12),Jtol=Jtol,
                                    max_iter=10000, verbose=False)
            unstable_param = [param/param_c for param in unstable_param]

            results[key]['unstable']['param'] = unstable_param
            results[key]['unstable']['infected_fraction'] = \
                    unstable_infected_fraction


    with open('./dat/Fig4_bif.pk', 'wb') as outfile:
        pickle.dump(results,outfile)
