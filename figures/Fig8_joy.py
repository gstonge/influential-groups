import numpy as np
import matplotlib.pyplot as plt
from gcm import *
import pickle

#structure
mmax = 4
gm = np.zeros(mmax+1)
gm[mmax] += 1

nmin = 2
nu = 1.5
gamma = 3.
nmax = 50
n= np.arange(nmax+1)
pn = np.zeros(nmax+1)
pn[nmin:] = n[nmin:]**(-gamma)
pn /= np.sum(pn)
state_meta = get_state_meta(mmax, nmax, gm, pn)

#infection
beta = lambda n,i,trate,nu: trate*i**nu
param_c = invasion_threshold_safe(beta, gm, pn,fixed_args=(nu,), max_param=0.1)
trate = param_c*1.1
inf_mat = infection_matrix(beta, nmax, args=(trate,nu))

#get stationary solution
r = stationary_state_safe(inf_mat, state_meta, r_upp=None, r_low=10**(-12))
v = state_from_mf(r,inf_mat,state_meta)
fni = v[mmax+1:].reshape((nmax+1,nmax+1))

with open('./dat/Fig6_joy.pk', 'wb') as outfile:
    pickle.dump(fni,outfile)
