import pickle
from gcm import *
import matplotlib.pyplot as plt

compute = True

#structure
mmin = 3
nmax = 5
mmax = 3
pn = np.zeros(nmax+1)
pn[4] += 1
gm = np.zeros(mmax+1)
gm[3] += 1
# state_meta = StateMeta(mmin, mmax, nmax, gm, pn)
state_meta = get_state_meta(mmax, nmax, gm, pn)
#infection
nu_list = [1.5, 1.7, 1.9080447119702004, 2.1, 2.3]
nu_c = 1.9080447119702004
beta = lambda n,i,trate,nu: trate*i**nu

if compute:
    #prepare result dict
    results = {nu:dict() for nu in nu_list}
    for nu in nu_list:
        results[nu]['stable'] = dict()
        if nu > nu_c:
            results[nu]['unstable'] = dict()
    #compute stable and unstable solutions
    for nu in nu_list:
        print(f'nu {nu}')
        param_c = invasion_threshold(beta,gm,pn,fixed_args=(nu,))
        param_init = 1.1*param_c
        #stable branch
        param_var = -0.005 #backward direction
        Jtol = 0.0001
        stable_param,stable_fixed_point,stable_infected_fraction = \
                stable_branch(beta,state_meta,param_init,param_var,
                              rtol=10**(-12), max_iter=20000,
                              fixed_args=(nu,),Jtol=Jtol,verbose=False)

        results[nu]['param_c'] = param_c
        results[nu]['stable']['param_list'] = list(stable_param/param_c)
        results[nu]['stable']['infected_fraction_list'] = list(
            stable_infected_fraction)
        if nu > nu_c:
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
            results[nu]['unstable']['param_list'] = list(
                unstable_param/param_c)
            results[nu]['unstable']['infected_fraction_list'] = list(
                unstable_infected_fraction)
    with open('./dat/Fig3_bif1.pk', 'wb') as outfile:
        pickle.dump(results,outfile)
