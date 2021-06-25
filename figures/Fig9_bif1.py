import pickle
from gcm import *
import matplotlib.pyplot as plt


#infection
nu = 0.5
beta = lambda n,i,trate,nu: trate*i**nu

#structure
# gamma_m = 3
# mmax = 100
# gm = np.zeros(mmax+1)
# gm[2:] += (np.arange(2,mmax+1)*1.)**(-gamma_m)
mmax = 4
gm = np.zeros(mmax+1)
gm[mmax] += 1
nmax = 50
gamma_n = 3
pn = np.zeros(nmax+1)
pn[2:] = (np.arange(2,nmax+1)*1.)**(-gamma_n)
pn /= np.sum(pn)
state_meta = get_state_meta(mmax, nmax, gm, pn)
nu_c = bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                              max_params=(1,1.8))
param_c = invasion_threshold(beta,gm,pn,fixed_args=(nu,))
print(nu_c,param_c)


#prepare result dict
results = dict()

#compute stable solutions
param_init = 5*param_c
param_var = -0.01*param_c #backward direction
Jtol = 0.0001
stable_param,stable_fixed_point,stable_infected_fraction = \
        stable_branch(beta,state_meta,param_init,param_var,
                      rtol=10**(-10), max_iter=3000,
                      fixed_args=(nu,),Jtol=Jtol,verbose=True)

#format fni to In per group
@jit(nopython=True)
def get_In_list(fni_list,nmax):
    ivec = np.arange(nmax+1,dtype=np.float64)
    In_list = np.zeros((len(fni_list),nmax+1))
    for j in range(len(fni_list)):
        for n in range(nmax+1):
            In_list[j][n] = np.sum(fni_list[j][n]*ivec/n)
    return In_list

fni_list = np.array([fni for _,fni in stable_fixed_point])
sm_list = np.array([sm for sm,_ in stable_fixed_point])
In_list = get_In_list(fni_list,nmax)

results['param_c'] = param_c
results['param_list'] = np.array(stable_param)
# results['fni_list'] = np.array(fni_list)
results['I_list'] = np.array(stable_infected_fraction)
results['In_list'] = np.array(In_list)

for n in [2,10,20,30,50]:
    plt.plot(results['param_list'],results['In_list'][:,n])
plt.axvline(param_c,0,1,color='k',ls='--')
plt.show()

with open('./dat/FigX_bif1.pk', 'wb') as outfile:
    pickle.dump(results,outfile)
