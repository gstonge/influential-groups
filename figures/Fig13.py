import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import loggamma
from gcm import *
from matplotlib import cm

compute = True

#infection
def beta(n,i,trate,epsilon):
    total = 0.
    if i > 0:
        total += epsilon
    if (i == n-1):
        total += trate-epsilon
    return total

nmax = 3
mmax = 6
gm = np.zeros(mmax+1)
pn = np.zeros(nmax+1)
pn[nmax] += 1
gm[mmax] += 1
state_meta = get_state_meta(mmax,nmax,gm,pn)

#plot parameters
font_size=8
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font',family='serif',serif='Computer Modern')
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)
plt.rc('axes', labelsize=font_size)
colorlist = ['#d0d1e6','#74a9cf','#045a8d']
width = 7.057/2
# height = width/(2*1.35)
height = width/2.

fig, ax = plt.subplots(1,1, figsize=(width, height))
plt.subplots_adjust(left=0.25, bottom=0.21, right=0.85, top=0.99,
                    wspace=0.20, hspace=0.25)


epsilon = [0.001,0.01,0.03,0.05]
color_list = ["#7fcdbb","#41b6c4","#1d91c0","#225ea8", "#0c2c84"]
param_range = [0.5,2.19]
for i,epsilon in enumerate(epsilon):
    #get stable branch
    trate_init = param_range[1]
    param_var = -0.05
    fixed_args = (epsilon,)
    param_list,state_list,Ilist = stable_branch(beta,state_meta,trate_init,
                                                param_var,fixed_args=fixed_args,
                                                verbose=False)
    ax.plot(param_list,Ilist,color=color_list[i],lw=1.5,
             label=fr"$\varepsilon = {epsilon}$")
    #get unstable branch
    fni = state_list[-1][1]
    param_init = param_list[-1]
    param_var = abs(param_list[-2]-param_list[-1])
    param_list,state_list,Ilist = unstable_branch(fni,beta,state_meta,param_init,
                                                  param_var,fixed_args=fixed_args,
                                                  verbose=False)
    ax.plot(param_list,Ilist,'--', color=color_list[i],lw=1.5)
    ax.plot([param_range[0],param_list[-1]],[0,0],color=color_list[i],lw=1.5)
plt.xlabel(r"$\lambda$")
plt.ylabel(r"Infected fraction $I^*$")
plt.legend(frameon=False)
plt.xlim(param_range)

plt.show()

