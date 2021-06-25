import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import brentq
from scipy.special import lambertw
from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size


#plot parameters
font_size=8
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font',family='serif',serif='Computer Modern')
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)
plt.rc('axes', labelsize=font_size)
plt.rc('legend', fontsize=font_size)


#color list
# color_list = ["#ef5310","#f59770","#fbdcd0","#a6dcf2","#4eb9e4","#1b86b1"]
# color_list = ["#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8", "#0c2c84"]
# color_list = ['#b35806','#f1a340','#fee0b6','#d8daeb','#998ec3','#542788']
color_list = ['#b2182b','#ef8a62','#fddbc7','#d1e5f0','#67a9cf','#2166ac']
newcm = LinearSegmentedColormap.from_list('ColorMap',
                                          list(reversed(color_list)))
blackcm = LinearSegmentedColormap.from_list('BlackColorMap',['#1a1a1a','#1a1a1a'])

#plot
width = 7.057
height = width/2
fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(width,height))

plt.subplots_adjust(left=0.085, bottom=0.12, right=0.88,
                    top=0.93, wspace=0.4, hspace=0.42)

#top row
#-------
filename_list = ['./dat/Fig10_hom_hom.pk','./dat/Fig10_hom_het.pk',
                 './dat/Fig10_het_het.pk']
#get the vmin and vmax
vmin = np.inf
vmax = -np.inf
for filename in filename_list:
    with open(filename, 'rb') as data:
        results = pickle.load(data)
    vmin = min([vmin,np.min(np.log10(results['zeta_arr']))])
    vmax = max([vmin,np.max(np.log10(results['zeta_arr']))])
vmax = max((abs(vmin),abs(vmax)))
vmin = -vmax

for ax,filename in zip(axes[0,:],filename_list):
    with open(filename, 'rb') as data:
        results = pickle.load(data)

    #find nu critical
    index = np.argmin(np.abs(np.array(results['zeta_lim'])-1.))
    ax.plot(results['nu_list'][index],np.amin(results['eps_arr']),'*',
            markersize=10,color='#1a1a1a',clip_on=False,zorder=10)

    CS = ax.contourf(results['nu_arr'],results['eps_arr'],
                     np.log10(results['zeta_arr']), cmap=newcm, zorder=-1,
                     vmin=vmin, vmax=vmax, levels=12)
    ax.contour(results['nu_arr'],results['eps_arr'],
               np.log10(results['zeta_arr']),levels=[0.],
               linestyles='dashed', cmap=blackcm, zorder=1)
    # ax.set_xlim([0.9,4.1])
    ax.set_xlabel(r"Nonlinear exponent $\nu$")
    ax.set_yscale('log')

pos = (0.0, 1.1) #for panel label

#left plot
#======================

axes[0,0].set_ylabel(r"Initial infected fraction $\epsilon$")
axes[0,0].text(pos[0], pos[1],
             r"(a) Homogeneous $g_m$ and $p_n$",
             fontsize=font_size,transform=axes[0,0].transAxes)


#center plot
#======================
axes[0,1].text(pos[0], pos[1],
             r"(b) Het. $g_m$ and hom. $p_n$",
             fontsize=font_size,transform=axes[0,1].transAxes)

#right plot
#======================
axes[0,2].text(pos[0], pos[1],
             r"(c) Heterogeneous $g_m$ and $p_n$",
             fontsize=font_size,transform=axes[0,2].transAxes)


x,y,w,h = axes[0,2].get_position().bounds
cax = fig.add_axes([0.90, y, 0.013, h])
cbar = fig.colorbar(CS, cax=cax)
cbar.ax.set_ylabel(r"$\log_{10} \zeta$")

#bottom row
#----------

strategies_label = {'spreaders':'Spreaders',
                   'groups':'Groups',
                   'random':'Random'}
strategies_color = {'spreaders':"#67a9cf",
                   'groups':"#ef8a62",
                   'random': "#666666"}
with open('./dat/Fig10_temp_evo.pk','rb') as filename:
    results = pickle.load(filename)
nu_list = list(sorted([nu for nu in results]))

panel_label_list = ['(d)','(e)','(f)']
for nu,ax,panel_label in zip(nu_list,axes[1],panel_label_list):
    t = results[nu]['t']
    for strategy,label in strategies_label.items():
        ax.plot(t, results[nu][strategy], label=label,
                color=strategies_color[strategy])
    if nu == 1:
        leg = ax.legend()
        leg.get_frame().set_linewidth(0.0)
        ax.set_ylabel(r"Infected fraction $I(t)$")
    ax.set_xlabel(r"Time $t$")
    ax.text(0.08, 0.8,panel_label,
             fontsize=font_size,transform=ax.transAxes,
             bbox=dict(facecolor='white',linewidth=0.,
                      alpha=0.9))

# marker_color = "#41b6c4"
marker_color = "#1a1a1a"
#put markers
axes[0,0].plot(1, 0.01, 'o', color=marker_color,mfc='white',mew=1)
axes[1,0].plot(0.25, 0.84, 'o', color=marker_color,mfc='white',mew=1,
         transform=axes[1,0].transAxes,zorder=10)
axes[0,0].plot(1.88, 0.01, 's', color=marker_color,mfc='white',mew=1)
axes[1,1].plot(0.25, 0.84, 's', color=marker_color,
         transform=axes[1,1].transAxes,mfc='white',mew=1,zorder=10)
axes[0,0].plot(3, 0.01, 'D', color=marker_color,mfc='white',mew=1)
axes[1,2].plot(0.25, 0.84, 'D', color=marker_color,
         transform=axes[1,2].transAxes,mfc='white',mew=1,zorder=10)

plt.savefig('../../paper/figs/Fig11.pdf')
plt.show()
