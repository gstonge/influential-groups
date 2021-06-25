import pickle
import matplotlib.pyplot as plt
import numpy as np
from labellines import labelLine, labelLines
from joyplot import *

#plot parameters
font_size=8
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font',family='serif',serif='Computer Modern')
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)
plt.rc('axes', labelsize=font_size)
plt.rc('legend', fontsize=font_size)
colorlist = ["#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8", "#0c2c84"]


#define figure
width = 7.057
height = width/4.2
nb_dist = 5
fig, axs = plt.subplots(ncols=3, nrows=nb_dist, figsize=(width,height))

#group the rows of two first columns
axes = []
for i in range(1):
    for j in range(2):
        gs = axs[i*nb_dist, j].get_gridspec()
        for ax in axs[i*nb_dist:(i+1)*nb_dist, j]:
            ax.remove()
        axbig = fig.add_subplot(gs[i*nb_dist:(i+1)*nb_dist, j])
        axes.append(axbig)

small_axes = axs[0:nb_dist,2]

plt.subplots_adjust(left=0.08, bottom=0.23, right=0.98, top=0.97, wspace=0.37,
                    hspace=0.25)


#Left side figure
#=====================
with open('./dat/Fig6_exp.pk', 'rb') as datafile:
    results = pickle.load(datafile)

first_element = next(results.__iter__())
gamma_list = list(results[first_element].keys())
nu_list = [nu for nu in results]

for i,nu in enumerate(nu_list):
    inv_threshold = [results[nu][gamma] for gamma in gamma_list]
    axes[0].plot(gamma_list,inv_threshold,'-', lw=1.5,
                     color=colorlist[2*i],
                 label=r"$\nu$="+f"{nu}")

leg = axes[0].legend(frameon=False)
# leg.get_frame().set_linewidth(0.0)
axes[0].set_xlabel(r"Group size exponent $\gamma_n$")
axes[0].set_ylabel(r"Invasion threshold $\lambda_\mathrm{c}$")
axes[0].set_ylim([0.001,0.165])
axes[0].text(0.9, 0.85, r'(a)', ha='center',
             va='center', transform=axes[0].transAxes,
             fontsize=font_size)


#middle figure
#=====================
with open('./dat/Fig6_size.pk', 'rb') as datafile:
    results = pickle.load(datafile)

first_element = next(results.__iter__())
nmax_list = list(results[first_element].keys())
nu_list = [nu for nu in results]

for i,nu in enumerate(nu_list):
    inv_threshold = [results[nu][nmax] for nmax in nmax_list]
    axes[1].loglog(nmax_list,inv_threshold,'-', lw=1.5,
                   color=colorlist[2*i],
                 label=r"$\nu$="+f"{nu}")
axes[1].scatter(50,results[1.5][50], 80, marker='*',color='#1a1a1a',zorder=3,
                         clip_on=False)

axes[1].loglog(nmax_list[20:],3*np.array(nmax_list[20:])**(-2.),'--', lw=1.5,
                   color='#282828',
                 label=r"$\sim n_\mathrm{max}^{-2}$")
leg = axes[1].legend(labelspacing=0.3, frameon=False)
# leg.get_frame().set_linewidth(0.0)
axes[1].set_xlabel(r"Maximal group size $n_\mathrm{max}$")
axes[1].set_ylabel(r"Invasion threshold $\lambda_\mathrm{c}$")
axes[1].text(0.9, 0.85, r'(b)', ha='center',
             va='center', transform=axes[1].transAxes,
             fontsize=font_size)

#Right side figure
#=====================
with open('./dat/Fig6_joy.pk', 'rb') as datafile:
    fni = pickle.load(datafile)

nlist = [10,20,30,40,50]
fni_dict = {fr"$n = {n}$":{i/n: fni[n][i] for i in range(n+1)} for n in nlist}
cm = LinearSegmentedColormap.from_list('cmap', ["#7fcdbb","#0c2c84"])
color_map = [cm(x) for x in np.linspace(0.8,0.1,nb_dist)]

joyplot(small_axes,fni_dict, xmin=-0.2, xmax=0.85, smooth=False,rescale=True,
        xticks=[0,0.2,0.4,0.6,0.8], color_map=color_map)

for ax in small_axes:
    ax.set_xticks([])
    l, b, w, h = ax.get_position().bounds
    ax.set_position([l, 0.09+0.65*b, w-0.01, h])
small_axes[-1].set_xlabel(r'Infected fraction $i/n$')
small_axes[-1].set_xticks([0,0.2,0.4,0.6,0.8])
small_axes[0].text(0., 1.65, r'(c)', fontsize=font_size,
                      transform=small_axes[0].transAxes)
small_axes[0].text(0.25,1.60,r'Distribution $f_{n,i}$',fontsize=font_size,
                      transform=small_axes[0].transAxes)
small_axes[0].scatter(0.13,1.8, 80, marker='*',color='#1a1a1a',zorder=3,
                     clip_on=False)

plt.savefig('../../paper/figs/Fig8.pdf')
# plt.show()


