import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import loggamma
from gcm import *


#Poisson membership and group size
nmax = 20
mmax = 20
m = np.arange(mmax+1)
n = np.arange(nmax+1)
pn = np.zeros(nmax+1)
param = 5
gm = np.exp(m*np.log(param) - loggamma(m+1))
gm[0:1] = 0
gm /= np.sum(gm)
pn = np.exp(n*np.log(param) - loggamma(n+1))
pn[0:2] = 0
pn /= np.sum(pn)
state_meta = get_state_meta(mmax, nmax, gm, pn)

#infection
beta = lambda n,i,trate,nu: trate*i**nu
nu_c = bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                              max_params=(1,10))

results = dict()
nu_list = [1,1.88,3]
tmax_list = [5,5,100]
for nu,tmax in zip(nu_list,tmax_list):
    #parameters
    trate_c = invasion_threshold_safe(beta,gm,pn, fixed_args=(nu,))
    trate = 3*trate_c
    epsilon = 10**(-2)
    inf_mat = infection_matrix(beta, nmax, args=(trate,nu))

    sm_S,fni_S = optimize_sm(epsilon,state_meta)
    sm_F,fni_F,_ = optimize_fni(epsilon, inf_mat, state_meta)
    sm_R,fni_R = initialize(state_meta, initial_density=epsilon)

    init_map = {'spreaders':(sm_S,fni_S),
                'groups':(sm_F,fni_F),
                'random':(sm_R,fni_R)}

    t = np.linspace(0,tmax,1000)
    results[nu] = dict()
    results[nu]['t'] = t
    for strategy in ['spreaders','groups','random']:
        sm,fni = init_map[strategy]
        v = flatten(sm,fni,state_meta)
        vvec = odeint(vector_field,v,t,args=(inf_mat,state_meta))
        Ivec = []
        for v in vvec:
            sm,fni = unflatten(v,state_meta)
            Ivec.append(1-np.sum(sm*gm))
        results[nu][strategy] = Ivec
        plt.plot(t,Ivec,label=strategy)
    plt.legend()
    plt.show()

#save data
with open('./dat/Fig10_temp_evo.pk', 'wb') as filename:
    pickle.dump(results,filename)
