
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

exp_list = ['coauth_exp2']
nu_list = [2]
case_list = ['sim_orig','sim_rew']
case_label_list = [r'\textbf{Original}',r'\textbf{Rewired}']
label_list = ['(b)','(c)']#,'(d)','(e)','(f)','(g)','(h)','(i)']

nb_exp = len(exp_list)
width = 7.057
height = nb_exp*width/3.8

fig, axes = plt.subplots(nb_exp,3, figsize=(width, height))
plt.subplots_adjust(left=0.08, bottom=0.22, right=0.98, top=0.88, wspace=0.30,
                    hspace=0.35)

for i,exp in enumerate(exp_list):
    #results from ame
    with open(f'../simulation/dat/{exp}_ame.pk', 'rb') as datafile:
        ame_results = pickle.load(datafile)
    param_list_ame_upper = np.array(
        ame_results['param_list_upper'])/ame_results['scale_c']
    param_list_ame_lower = np.array(
        ame_results['param_list_lower'])/ame_results['scale_c']

    #structure data
    datafile = '../simulation/socio_data/form_coauth-dblp_sub.pk'
    with open(datafile,'rb') as filename:
        data = pickle.load(filename)
    gm = data['gm']
    pn = data['pn']
    mmax = data['mmax']
    nmax = data['nmax']

    #plot structure data
    if nb_exp > 1:
        ax = axes[i][0]
    else:
        ax = axes[0]
    ax.text(0.5, 1.1, r'\textbf{Structure}', ha='center',va='center',
                transform=ax.transAxes,fontsize=font_size)
    axins = inset_axes(ax, width="45%", height="45%",borderpad=0)

    ax.scatter(np.arange(mmax+1)[1:],gm[1:], 18, color='#b2b2b2')
    pl_x = np.arange(6,150)
    pl_y = 1.3*pl_x**(-2.3)
    ax.plot(pl_x,pl_y, color='#1a1a1a', label=r"$m^{-2.3}$")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([6*10**-6,20])
    ax.set_xlim([0.5,10**4])
    ax.set_ylabel("$g_m$")
    ax.set_xlabel("$m$")
    ax.text(0.12, 0.85, '(a)', ha='center',va='center',
            transform=ax.transAxes,fontsize=font_size)
    ax.legend(frameon=False,loc='lower left')


    axins.bar(np.arange(nmax+1)[2:],pn[2:], facecolor="#41b6c4",
             edgecolor="#41b6c4",width=1.)
    axins.set_ylabel("$p_n$")
    axins.set_xlabel("$n$")
    axins.set_xticks([10,20])
    axins.set_yscale('log')

    #plot simulation vs ame
    for j,case in enumerate(case_list):
        #simulation data
        with open(f'../simulation/dat/{exp}_{case}.pk', 'rb') as datafile:
            sim_results = pickle.load(datafile)
        param_list_sim = np.array(sim_results['param_list'])/\
                ame_results['scale_c']
        if nb_exp > 1:
            ax = axes[i][j+1]
        else:
            ax = axes[j+1]
        ax.plot(param_list_ame_upper, ame_results['I_list_upper'],
                color="#0c2c84")
        ax.plot(param_list_ame_lower, ame_results['I_list_lower'],
                color="#0c2c84")
        ax.errorbar(param_list_sim,sim_results['mean_I_list_upper'],
                     yerr=sim_results['std_I_list_upper'],
                     color="#7fcdbb", fmt='o', ms=5,elinewidth=1.)
        if i == nb_exp-1:
            ax.set_xlabel(r"$\lambda/\lambda_\mathrm{c}$")
        # if j == 0:
        ax.text(0.3, 0.85, fr'$\nu = {nu_list[i]}$', ha='center',va='center',
                transform=ax.transAxes,fontsize=font_size)
        ax.set_ylabel(r"Infected fraction $I^*$")

        ax.text(0.12, 0.85, f'{label_list[3*i+j]}', ha='center',va='center',
                transform=ax.transAxes,fontsize=font_size)
        if i == 0:
            ax.text(0.5, 1.1, case_label_list[j], ha='center',va='center',
                    transform=ax.transAxes,fontsize=font_size)

plt.show()
