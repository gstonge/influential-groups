import pickle
import matplotlib.pyplot as plt
import numpy as np
from labellines import labelLine, labelLines
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#plot parameters
font_size=8
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font',family='serif',serif='Computer Modern')
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)
plt.rc('axes', labelsize=font_size)
colorlist = ["#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8", "#0c2c84"]

nu = 2

width = 7.057/2
height = width/1.8

fig, ax = plt.subplots(1,1, figsize=(width, height))
plt.subplots_adjust(left=0.22, bottom=0.2, right=0.90, top=0.98, wspace=0.30,
                    hspace=0.35)

#results from ame
with open(f'../simulation/dat/modified_coauth_exp2_ame.pk', 'rb') as datafile:
    ame_results = pickle.load(datafile)
param_list_ame_upper = np.array(
    ame_results['param_list_upper'])/ame_results['scale_c']
param_list_ame_lower = np.array(
    ame_results['param_list_lower'])/ame_results['scale_c']

#structure data
datafile = '../simulation/socio_data/form_modified_coauth-dblp_sub.pk'
with open(datafile,'rb') as filename:
    data = pickle.load(filename)
gm = data['gm']
pn = data['pn']
mmax = data['mmax']
nmax = data['nmax']


#simulation data
with open(f'../simulation/dat/modified_coauth_exp2_sim_rew.pk', 'rb') as datafile:
    sim_results = pickle.load(datafile)
param_list_sim = np.array(sim_results['param_list'])/\
        ame_results['scale_c']
ax.plot(param_list_ame_upper, ame_results['I_list_upper'],
        color="#0c2c84")
ax.plot(param_list_ame_lower, ame_results['I_list_lower'],
        color="#0c2c84")
ax.errorbar(param_list_sim,sim_results['mean_I_list_upper'],
             yerr=sim_results['std_I_list_upper'],
             color="#7fcdbb", fmt='o', ms=5,elinewidth=1.)
ax.set_xlabel(r"$\lambda/\lambda_\mathrm{c}$")
ax.text(0.2, 0.85, fr'$\nu = {nu}$', ha='center',va='center',
        transform=ax.transAxes,fontsize=font_size)
ax.set_ylabel(r"Infected fraction $I^*$")


plt.show()
