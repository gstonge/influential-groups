import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import loggamma
from gcm import *
from matplotlib.colors import LinearSegmentedColormap

color_list = ["#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8", "#0c2c84"]
newcm = LinearSegmentedColormap.from_list('ColorMap',
                                          list(reversed(color_list[:-1])))
blackcm = LinearSegmentedColormap.from_list('BlackColorMap',['#1a1a1a','#1a1a1a'])

#Poisson membership and group size
nmax = 20
mmax = 100
m = np.arange(mmax+1)
n = np.arange(nmax+1)
gm = np.zeros(mmax+1)
gm[2:] = (m[2:]*1.)**(-3.)
gm /= np.sum(gm)
param = 5
pn = np.exp(n*np.log(param) - loggamma(n+1))
pn[0:2] = 0
pn /= np.sum(pn)
state_meta = get_state_meta(mmax, nmax, gm, pn)

#infection
trate = 1 #independent from this
beta = lambda n,i,trate,nu: trate*i**nu
nu_c = bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                              max_params=(3,8))

#parameters
nu_list = np.linspace(0.9,4.1,100)
epsilon_list = np.logspace(-5,-1,100)

#get limit zeta
zeta_lim = []
for nu in nu_list:
    inf_mat = infection_matrix(beta, nmax, args=(trate,nu))
    zeta_lim.append(zeta_limit(inf_mat,state_meta))
    print(nu,zeta_lim[-1])

#get zeta arr
zeta_arr = np.zeros((len(epsilon_list),len(nu_list)))
zeta_arr[:] = np.nan
for i,epsilon in enumerate(epsilon_list):
    for j,nu in enumerate(nu_list):
        inf_mat = infection_matrix(beta, nmax, args=(trate,nu))
        sm_S,fni_S = optimize_sm(epsilon,state_meta)
        sm_F,fni_F,_ = optimize_fni(epsilon, inf_mat, state_meta)
        Phi_F = objective_function(fni_F,inf_mat,state_meta)
        Phi_S = objective_function(fni_S,inf_mat,state_meta)
        if Phi_F > 0 and Phi_S > 0:
            zeta_arr[i,j] = Phi_F/Phi_S

nu_arr, eps_arr = np.meshgrid(nu_list, epsilon_list)
plt.contourf(nu_arr, eps_arr, np.log10(zeta_arr), cmap=newcm, zorder=-1,
            levels=13)
plt.yscale('log')
plt.colorbar()
plt.contour(nu_arr, eps_arr, np.log10(zeta_arr), cmap=blackcm, levels=[0.])
plt.show()

#save data
results = dict()
results['nu_arr'] = nu_arr
results['eps_arr'] = eps_arr
results['zeta_arr'] = zeta_arr
results['nu_list'] = nu_list
results['zeta_lim'] = zeta_lim

with open('./dat/Fig10_hom_het.pk', 'wb') as filename:
    pickle.dump(results,filename)
