import pickle
import matplotlib.pyplot as plt
from labellines import labelLine, labelLines
from gcm import *
import numpy as np

#plot parameters
font_size=8
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font',family='serif',serif='Computer Modern')
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)
plt.rc('axes', labelsize=font_size)
colorlist = list(reversed(["#c7e9b4","#41b6c4","#1d91c0","#0c2c84"]))
width = 7.057
height = width/3.7

fig, axes = plt.subplots(1,3, figsize=(width, height))
plt.subplots_adjust(left=0.07, bottom=0.22, right=0.98, top=0.88, wspace=0.30,
                    hspace=0.25)

nlist = list(reversed([2,10,30,50]))
color_map = {n:color for n,color in zip(nlist,colorlist)}

#Left side figure
#=====================
with open('./dat/FigX_bif1.pk', 'rb') as datafile:
    results = pickle.load(datafile)

for n in nlist:
    axes[0].plot([0.6,1],[0,0],
             '-', lw=1.5, color=color_map[n])
    axes[0].plot(results['param_list']/results['param_c'],
                 results['In_list'][:,n],
             '-', lw=1.5, color=color_map[n],
             label=rf"$n = {n}$")
# axes[1].set_xticks([0.95,1.,1.05,1.1])
# axes[1].set_yticks([0.,0.2,0.4,0.6])
axes[0].set_xlabel(r"$\lambda/\lambda_\mathrm{c}$")
axes[0].set_ylabel(r"Infected fraction $I_n^*$")
xvals = [1.,1.03,1.06]
# labelLines(list(axes[1].get_lines()), zorder=2.5,align=True,xvals=xvals,
           # color='black',fontsize=font_size)
axes[0].text(0.12, 1.1, r'(a) $\nu = 0.5,\; g_m = \delta_{m,4}$', ha='left',
             va='center', transform=axes[0].transAxes,
            fontsize=font_size)
leg = axes[0].legend()
leg.get_frame().set_linewidth(0.0)


#middle side figure
#=====================

with open('./dat/FigX_bif2.pk', 'rb') as datafile:
    results = pickle.load(datafile)

for n in nlist:
    axes[1].plot([0.6,1],[0,0],
             '-', lw=1.5, color=color_map[n])
    axes[1].plot(results['param_list']/results['param_c'],
                 results['In_list'][:,n],
             '-', lw=1.5, color=color_map[n],
             label=rf"$n = {n}$")
# axes[1].set_xticks([0.95,1.,1.05,1.1])
axes[1].set_yticks([0.,0.2,0.4,0.6,0.8])
axes[1].set_xlabel(r"$\lambda/\lambda_\mathrm{c}$")
axes[1].set_ylabel(r"Infected fraction $I_n^*$")
xvals = [1.,1.03,1.06]
# labelLines(list(axes[1].get_lines()), zorder=2.5,align=True,xvals=xvals,
           # color='black',fontsize=font_size)
axes[1].text(0.12, 1.1, r'(b) $\nu = 1.5,\; g_m = \delta_{m,4}$', ha='left',
             va='center', transform=axes[1].transAxes,
            fontsize=font_size)


#right panel
#-----------

with open('./dat/FigX_bif3.pk', 'rb') as datafile:
    results = pickle.load(datafile)

for n in nlist:
    axes[2].plot([0.95,1],[0,0],
             '-', lw=1.5, color=color_map[n],zorder=-10,alpha=0.2)
    axes[2].plot(results['stable']['param_list']/results['param_c'],
                 results['stable']['In_list'][:,n],
             '-', lw=1.5, color=color_map[n],
             label=rf"$n = {n}$", alpha = 0.2)
    axes[2].plot(results['unstable']['param_list']/results['param_c'],
             results['unstable']['In_list'][:,n],
             '--', lw=1.5, color=color_map[n])
# axes[1].set_xticks([0.95,1.,1.05,1.1])
# axes[1].set_yticks([0.,0.2,0.4,0.6])
# axes[2].set_xlim([0.94,1.03])
# axes[2].set_ylim([-0.04,0.72])
axes[2].set_xlabel(r"$\lambda/\lambda_\mathrm{c}$")
axes[2].set_ylabel(r"Infected fraction $I_n^*$")
xvals = [1.,1.03,1.06]
# labelLines(list(axes[2].get_lines()), zorder=2.5,align=True,xvals=xvals,
           # color='black',fontsize=font_size)
axes[2].text(0.12, 1.1, r'(c) $\nu = 1.5,\; g_m = \delta_{m,10}$', ha='left',
             va='center', transform=axes[2].transAxes,
            fontsize=font_size)



plt.savefig('../../paper/figs/Fig9.pdf')
plt.show()


