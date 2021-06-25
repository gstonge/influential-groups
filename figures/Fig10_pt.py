
import pickle
from gcm import *

compute = True

#infection
beta = lambda n,i,trate,nu: trate*i**nu

#structures
mmax = 4
gm = np.zeros(mmax+1)
gm[mmax] += 1

nmin = 2
nmax = 15
n = np.arange(nmax+1)
gamma = 3.
pn = np.zeros(nmax+1)
pn[nmin:] = n[nmin:]**(-gamma)
pn[0:nmin] = 0.
pn /= np.sum(pn)

state_meta_1 = get_state_meta(mmax, nmax, gm, pn)


nmax = 40
n = np.arange(nmax+1)
pn = np.zeros(nmax+1)
pn[nmin:] = n[nmin:]**(-gamma)
pn[0:nmin] = 0.
pn /= np.sum(pn)
state_meta_2 = get_state_meta(mmax, nmax, gm, pn)

state_meta_dict = {"struct1":state_meta_1,
                   "struct2":state_meta_2}


if compute:
    #prepare result dict
    results = {key:dict() for key in state_meta_dict}
    #compute solution
    for key,state_meta in state_meta_dict.items():
        gm = state_meta[3]
        pn = state_meta[4]

        nu_c = bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                                      max_params=(1,3))
        param_c = invasion_threshold(beta,gm,pn,fixed_args=(nu_c,))
        print(key)
        print(nu_c,param_c)

        param_init = 0.2
        #stable branch
        param_var = -0.01*param_c #backward direction
        Jtol = 0.0001
        stable_param,stable_fixed_point,stable_infected_fraction = \
                stable_branch(beta,state_meta,param_init,param_var,
                              rtol=10**(-10), max_iter=3000,
                              fixed_args=(nu_c,),Jtol=Jtol,verbose=True)

        results[key]['param_list'] = list(stable_param)
        results[key]['infected_fraction_list'] = list(stable_infected_fraction)

    with open('./dat/Fig7_pt.pk', 'wb') as outfile:
        pickle.dump(results,outfile)
