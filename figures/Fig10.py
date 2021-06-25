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
width = 7.057
height = width/4

fig, axes = plt.subplots(1,3, figsize=(width, height))
plt.subplots_adjust(left=0.07, bottom=0.23, right=0.98, top=0.97, wspace=0.30,
                    hspace=0.25)


with open('./dat/Fig7_bist_thresh.pk', 'rb') as datafile:
    results = pickle.load(datafile)

#Left side figure
#=====================
first_element = next(results.__iter__())
gamma_list = list(results[first_element].keys())
nmax_list = [nmax for nmax in results]

for i,nmax in enumerate([20,40]):
    bthreshold = [results[nmax][gamma] for gamma in gamma_list]
    axes[0].plot(gamma_list,bthreshold,'-', lw=1.5,
                     color=colorlist[2*i+2],
                 label=r"$n_\mathrm{max}$="+f"{nmax}")

leg = axes[0].legend(loc='lower right')
leg.get_frame().set_linewidth(0.0)
axes[0].set_xlabel(r"Group exponent $\gamma_n$")
axes[0].set_ylabel(r"Bistability threshold $\nu_\mathrm{c}$")
# axes[0].set_xlim([1.8,5.35])

axes[0].text(0.12, 0.85, r'(a)', ha='center',
             va='center', transform=axes[0].transAxes,
             fontsize=font_size)

#middle figure
#=====================
plot_gamma = [2.1,3.,4.]
for i,gamma in enumerate(plot_gamma):
    bthreshold = [results[nmax][gamma] for nmax in nmax_list]
    axes[1].plot(nmax_list,bthreshold,'-', lw=1.5,
                     color=colorlist[2*i],
                 label=r"$\gamma_n$="+f"{gamma}")

leg = axes[1].legend(loc='upper right')
leg.get_frame().set_linewidth(0.0)
axes[1].set_xlabel(r"Maximal group size $n_\mathrm{max}$")
axes[1].set_ylabel(r"Bistability threshold $\nu_\mathrm{c}$")

axes[1].text(0.12, 0.85, r'(b)', ha='center',
             va='center', transform=axes[1].transAxes,
             fontsize=font_size)
axes[1].set_ylim([0.98,4.8])

#Right side figure
#=====================
with open('./dat/Fig7_pert.pk', 'rb') as datafile:
    results = pickle.load(datafile)

for key in reversed(["struct1","struct2"]):
    structure = results[key]
    if key == "struct1":
        label = r"Regular"
        color = colorlist[1]
    elif key == "struct2":
        label = r"Perturbed"
        color = colorlist[4]

    for stability,solution in structure.items():
        if stability == "stable":
            axes[2].plot(solution['param_list'],
                         solution['infected_fraction_list'],'-', lw=1.5,
                         color=color,label=label)
        elif stability == "unstable":
            axes[2].plot(solution['param_list'],
                         solution['infected_fraction_list'],'--', lw=1.5,
                         color=color)
            axes[2].plot((0,solution['param_list'][-1]),(0,0),'-', lw=1.5,
                         color=color)
# axes[2].set_xticks([0.95,1.,1.05,1.1])
axes[2].set_xlabel(r"$\lambda$")
axes[2].set_ylabel(r"Infected fraction $I^*$")
# xvals = [1.,1.06]
# labelLines(list(axes[2].get_lines()), zorder=2.5,align=True,#xvals=xvals,
           # color='black',fontsize=font_size)
axes[2].text(0.12, 0.85, r'(c)', ha='center',
             va='center', transform=axes[2].transAxes,
            fontsize=font_size)
axes[2].legend(frameon=False)


plt.savefig('../../paper/figs/Fig9.pdf')
# plt.show()


