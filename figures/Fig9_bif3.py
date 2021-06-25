import pickle
from gcm import *
import matplotlib.pyplot as plt


#infection
nu = 1.5
beta = lambda n,i,trate,nu: trate*i**nu

#structure
mmax = 10
gm = np.zeros(mmax+1)
gm[mmax] += 1
nmax = 50
gamma_n = 3
pn = np.zeros(nmax+1)
pn[2:] = (np.arange(2,nmax+1)*1.)**(-gamma_n)
pn /= np.sum(pn)
state_meta = get_state_meta(mmax, nmax, gm, pn)
# nu_c = bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                              # max_params=(1,7))
param_c = invasion_threshold(beta,gm,pn,fixed_args=(nu,))
print(param_c)


#prepare result dict
results = {'stable':dict(), 'unstable':dict()}

#compute stable solutions
param_init = 1.03*param_c
param_var = -0.01*param_c #backward direction
Jtol = 0.0001
stable_param,stable_fixed_point,stable_infected_fraction = \
        stable_branch(beta,state_meta,param_init,param_var,
                      rtol=10**(-10), max_iter=3000,
                      fixed_args=(nu,),Jtol=Jtol,verbose=True)

#get unstable solution
param_var = abs(stable_param[-2]-stable_param[-1])
fni = stable_fixed_point[-1][1]
param_init = stable_param[-1]
unstable_param,unstable_fixed_point,unstable_infected_fraction = \
        unstable_branch(fni,beta,state_meta,param_init,
                        param_var,fixed_args=(nu,),init_iter=100,
                        h=10**(-2),
                        rtol=10**(-12),Jtol=Jtol,
                        max_iter=10000, verbose=True)

#format fni to In per group
@jit(nopython=True)
def get_In_list(fni_list,nmax):
    ivec = np.arange(nmax+1,dtype=np.float64)
    In_list = np.zeros((len(fni_list),nmax+1))
    for j in range(len(fni_list)):
        for n in range(nmax+1):
            In_list[j][n] = np.sum(fni_list[j][n]*ivec/n)
    return In_list

fni_list_stable = np.array([fni for _,fni in stable_fixed_point])
sm_list_stable = np.array([sm for sm,_ in stable_fixed_point])
In_list_stable = get_In_list(fni_list_stable,nmax)

fni_list_unstable = np.array([fni for _,fni in unstable_fixed_point])
sm_list_unstable = np.array([sm for sm,_ in unstable_fixed_point])
In_list_unstable = get_In_list(fni_list_unstable,nmax)

results['param_c'] = param_c
results['stable']['param_list'] = np.array(stable_param)
# results['stable']['fni_list'] = np.array(fni_list_stable)
results['stable']['I_list'] = np.array(stable_infected_fraction)
results['stable']['In_list'] = np.array(In_list_stable)

results['unstable']['param_list'] = np.array(unstable_param)
# results['unstable']['fni_list'] = np.array(fni_list_unstable)
results['unstable']['I_list'] = np.array(unstable_infected_fraction)
results['unstable']['In_list'] = np.array(In_list_unstable)


for n in [2,10,20,30,50]:
    plt.plot(results['stable']['param_list'],results['stable']['In_list'][:,n])
    plt.plot(results['unstable']['param_list'],results['unstable']['In_list'][:,n])
plt.axvline(param_c,0,1,color='k',ls=':')
plt.show()

with open('./dat/FigX_bif3.pk', 'wb') as outfile:
    pickle.dump(results,outfile)
