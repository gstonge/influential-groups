import pickle
import matplotlib.pyplot as plt
import numpy as np
from labellines import labelLine, labelLines

#plot parameters
font_size=8
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font',family='serif',serif='Computer Modern')
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)
plt.rc('axes', labelsize=font_size)
colorlist = ["#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8", "#0c2c84"]

exp_list = ['Lyon_exp1','Lyon_exp4']#['Lyon_exp0.5','Lyon_exp1','Lyon_exp4']
nu_list = [1,4] #[0.5,1,4]
# exp_list = ['Thiers_exp1','Thiers_exp4']
case_list = ['sim_orig','sim_rew','sim_rew-exp']
case_label_list = [r'\textbf{Original}',r'\textbf{Rewired}',
                 r'\textbf{Rewired + expanded}']
label_list = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)']

nb_exp = len(exp_list)
width = 7.057
height = nb_exp*width/4.3

fig, axes = plt.subplots(nb_exp,3, figsize=(width, height))
plt.subplots_adjust(left=0.07, bottom=0.12, right=0.98, top=0.94, wspace=0.30,
                    hspace=0.35)

for i,exp in enumerate(exp_list):
    #results from ame
    with open(f'../simulation/dat/{exp}_ame.pk', 'rb') as datafile:
        ame_results = pickle.load(datafile)
    param_list_ame_upper = np.array(
        ame_results['param_list_upper'])/ame_results['scale_c']
    param_list_ame_lower = np.array(
        ame_results['param_list_lower'])/ame_results['scale_c']

    for j,case in enumerate(case_list):
        #simulation data
        with open(f'../simulation/dat/{exp}_{case}.pk', 'rb') as datafile:
            sim_results = pickle.load(datafile)
        param_list_sim = np.array(sim_results['param_list'])/\
                ame_results['scale_c']
        if nb_exp > 1:
            ax = axes[i][j]
        else:
            ax = axes[j]
        ax.plot(param_list_ame_upper, ame_results['I_list_upper'],
                color="#0c2c84")
        ax.plot(param_list_ame_lower, ame_results['I_list_lower'],
                color="#0c2c84")
        ax.errorbar(param_list_sim,sim_results['mean_I_list_upper'],
                     yerr=sim_results['std_I_list_upper'],
                     color="#7fcdbb", fmt='o', ms=5,elinewidth=1.)
        if i > 0:
            param_list_ame_unstable = np.array(
                ame_results['param_list_unstable'])/ame_results['scale_c']
            ax.plot(param_list_ame_unstable, ame_results['I_list_unstable'],
                    color="#0c2c84", ls='--')

        if j == 2 and i > 0:

            ax.errorbar(param_list_sim,sim_results['mean_I_list_lower'],
                         yerr=sim_results['std_I_list_lower'],
                         color="#1d91c0", fmt='^', ms=5,elinewidth=1.)
        if i == nb_exp-1:
            ax.set_xlabel(r"$\lambda/\lambda_\mathrm{c}$")
        if j == 0:
            ax.text(0.3, 0.85, fr'$\nu = {nu_list[i]}$', ha='center',va='center',
                transform=ax.transAxes,fontsize=font_size)
            ax.set_ylabel(r"Infected fraction $I^*$")

        ax.text(0.12, 0.85, f'{label_list[3*i+j]}', ha='center',va='center',
                transform=ax.transAxes,fontsize=font_size)
        if i == 0:
            ax.text(0.5, 1.1, case_label_list[j], ha='center',va='center',
                    transform=ax.transAxes,fontsize=font_size)

plt.show()
