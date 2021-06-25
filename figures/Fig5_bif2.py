import pickle
from gcm import *

compute = True

#structures
nmax = 3
mmax = 6
pn = np.zeros(nmax+1)
pn[nmax] += 1
gm = np.zeros(mmax+1)
gm[mmax] += 1
state_meta_1 = get_state_meta(mmax, nmax, gm, pn)

nmax = 4
mmax = 4
pn = np.zeros(nmax+1)
pn[nmax] += 1
gm = np.zeros(mmax+1)
gm[mmax] += 1
state_meta_2 = get_state_meta(mmax, nmax, gm, pn)


nmax = 5
mmax = 3
pn = np.zeros(nmax+1)
pn[nmax] += 1
gm = np.zeros(mmax+1)
gm[mmax] += 1
state_meta_3 = get_state_meta(mmax, nmax, gm, pn)

# #structures
# nmax = 4
# mmax = 4
# pn = np.zeros(nmax+1)
# pn[nmax] += 1
# gm = np.zeros(mmax+1)
# gm[mmax] += 1
# state_meta_1 = get_state_meta(mmax, nmax, gm, pn)

# nmax = 4
# mmax = 5
# pn = np.zeros(nmax+1)
# pn[nmax] += 1
# gm = np.zeros(mmax+1)
# gm[mmax] += 1
# state_meta_2 = get_state_meta(mmax, nmax, gm, pn)


# nmax = 4
# mmax = 6
# pn = np.zeros(nmax+1)
# pn[nmax] += 1
# gm = np.zeros(mmax+1)
# gm[mmax] += 1
# state_meta_3 = get_state_meta(mmax, nmax, gm, pn)


state_meta_dict = {"struct1":state_meta_1,
                   "struct2":state_meta_2,
                   "struct3":state_meta_3}

#infection
nu = 1.8
beta = lambda n,i,trate,nu: trate*i**nu

if compute:
    #prepare result dict
    results = {key:dict() for key in state_meta_dict}
    for key in state_meta_dict:
        results[key]['stable'] = dict()
        if key == "struct3":
            results[key]['unstable'] = dict()
    #compute stable and unstable solutions
    for key,state_meta in state_meta_dict.items():
        gm = state_meta[3]
        pn = state_meta[4]

        param_c = invasion_threshold(beta,gm,pn,fixed_args=(nu,))
        param_init = 1.1*param_c
        #stable branch
        param_var = -0.005 #backward direction
        Jtol = 0.0001
        stable_param,stable_fixed_point,stable_infected_fraction = \
                stable_branch(beta,state_meta,param_init,param_var,
                              rtol=10**(-12), max_iter=20000,
                              fixed_args=(nu,),Jtol=Jtol,verbose=False)

        results[key]['stable']['param_list'] = list(stable_param/param_c)
        results[key]['stable']['infected_fraction_list'] = list(
            stable_infected_fraction)
        if key == "struct3":
            #get unstable branch
            param_var = abs(stable_param[-2]-stable_param[-1])*param_c
            fni = stable_fixed_point[-1][1]
            param_init = stable_param[-1]
            unstable_param,unstable_fixed_point,unstable_infected_fraction = \
                    unstable_branch(fni,beta,state_meta,param_init,
                                    param_var,fixed_args=(nu,),init_iter=100,
                                    h=10**(-2),
                                    rtol=10**(-12),Jtol=Jtol,
                                    max_iter=10000, verbose=False)
            results[key]['unstable']['param_list'] = list(
                unstable_param/param_c)
            results[key]['unstable']['infected_fraction_list'] = list(
                unstable_infected_fraction)
    with open('./dat/Fig3_bif2.pk', 'wb') as outfile:
        pickle.dump(results,outfile)
