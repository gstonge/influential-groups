import pickle
from gcm import *

compute = True

#infection
nu = 2.3
beta = lambda n,i,trate,nu: trate*i**nu

#structures
nmax = 4
mmax = 4
pn = np.zeros(nmax+1)
pn[nmax] += 1
gm = np.zeros(mmax+1)
gm[mmax] += 1
state_meta_1 = get_state_meta(mmax, nmax, gm, pn)
print(bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                              max_params=(1,7)))
print(invasion_threshold(beta,gm,pn,fixed_args=(nu,)))


epsilon = 10**(-3)
nmin = 4
nmax = 15
mmax = 4
pn = np.zeros(nmax+1)
pn[nmin] += 1-epsilon
pn[nmax] += epsilon
gm = np.zeros(mmax+1)
gm[mmax] += 1
state_meta_2 = get_state_meta(mmax, nmax, gm, pn)
print(bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                              max_params=(1,7)))
print(invasion_threshold(beta,gm,pn,fixed_args=(nu,)))

state_meta_dict = {"struct1":state_meta_1,
                   "struct2":state_meta_2}



if compute:
    #prepare result dict
    results = {key:dict() for key in state_meta_dict}
    for key in state_meta_dict:
        results[key]['stable'] = dict()
        if key == "struct1":
            results[key]['unstable'] = dict()
    #compute stable and unstable solutions
    for key,state_meta in state_meta_dict.items():
        gm = state_meta[3]
        pn = state_meta[4]
        param_init = 0.2
        param_c = invasion_threshold(beta,gm,pn,fixed_args=(nu,))
        #stable branch
        param_var = -0.01*param_c #backward direction
        Jtol = 0.0001
        stable_param,stable_fixed_point,stable_infected_fraction = \
                stable_branch(beta,state_meta,param_init,param_var,
                              rtol=10**(-10), max_iter=3000,
                              fixed_args=(nu,),Jtol=Jtol,verbose=True)

        results[key]['stable']['param_list'] = list(stable_param)
        results[key]['stable']['infected_fraction_list'] = list(
            stable_infected_fraction)
        if key == "struct1":
            #get unstable branch
            param_var = abs(stable_param[-2]-stable_param[-1])*param_c
            fni = stable_fixed_point[-1][1]
            param_init = stable_param[-1]
            unstable_param,unstable_fixed_point,unstable_infected_fraction = \
                    unstable_branch(fni,beta,state_meta,param_init,
                                    param_var,fixed_args=(nu,),init_iter=100,
                                    h=10**(-2),
                                    rtol=10**(-12),Jtol=Jtol,
                                    max_iter=3000, verbose=True)
            results[key]['unstable']['param_list'] = list(
                unstable_param)
            results[key]['unstable']['infected_fraction_list'] = list(
                unstable_infected_fraction)
    with open('./dat/Fig7_pert.pk', 'wb') as outfile:
        pickle.dump(results,outfile)
