import pickle
import matplotlib.pyplot as plt
import numpy as np

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

fig, axes = plt.subplots(1,2, figsize=(width, height))
plt.subplots_adjust(left=0.17, bottom=0.22, right=0.88, top=0.98, wspace=0.45,
                    hspace=0.25)


#Left side figure
#=====================
with open('./dat/Fig4_bif.pk', 'rb') as datafile:
    results = pickle.load(datafile)

for key in reversed(["struct1","struct2"]):
    if key == "struct1":
        label = r"$\left < m^3 \right > \approx 206 $"
        color = colorlist[-2]
    elif key == "struct2":
        label = r"$\left < m^3 \right > \approx 267$"
        color = colorlist[1]
    structure = results[key]
    gm = structure['gm']
    cum_gm = np.cumsum(gm)
    mmin = np.searchsorted(cum_gm, 0, side='right')
    m = np.arange(len(gm))
    axes[0].loglog(m[mmin:],gm[mmin:],'-', lw=1.5,
                 color=color,label=label)


axes[0].set_ylim([2*10**(-8),20])
axes[0].set_xlabel(r"membership $m$")
axes[0].set_ylabel(r"Distribution $g_m$")
axes[0].legend(frameon=False)
axes[0].text(0.12, 0.85, r'(a)', ha='center',
             va='center', transform=axes[0].transAxes,
             fontsize=font_size)
#Right side figure
#=====================
with open('./dat/Fig4_bif.pk', 'rb') as datafile:
    results = pickle.load(datafile)

for key in reversed(["struct1","struct2"]):
    structure = results[key]
    if key == "struct1":
        label = r"$\left < m^3 \right > = $"
        color = colorlist[-2]
    elif key == "struct2":
        label = r"$\left < m^3 \right > =$"
        color = colorlist[1]

    for stability,solution in structure.items():
        if stability == "stable":
            axes[1].plot(solution['param'],
                         solution['infected_fraction'],'-', lw=1.5,
                         color=color)
        elif stability == "unstable":
            axes[1].plot(solution['param'],
                         solution['infected_fraction'],'--', lw=1.5,
                         color=color)
            axes[1].plot((0.96,solution['param'][-1]),(0,0),'-', lw=1.5,
                         color=color)

axes[1].set_xticks([0.95,1.,1.05,1.1])
axes[1].set_xlabel(r"$\lambda/\lambda_\mathrm{c}$")
axes[1].set_ylabel(r"Infected fraction $I^*$")
axes[1].text(0.12, 0.85, r'(b)', ha='center',
             va='center', transform=axes[1].transAxes,
             fontsize=font_size)
plt.savefig('../../paper/figs/Fig6.pdf')
plt.show()


