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
color_list = ["#67a9cf", "#ef8a62"]

width = 7.057
height = width/3.5

fig, axes = plt.subplots(1,2, figsize=(width, height))
plt.subplots_adjust(left=0.1, bottom=0.2, right=0.97, top=0.90, wspace=0.35,
                    hspace=0.45)

with open('../simulation/dat/Lyon_mixing.pk','rb') as filename:
    results = pickle.load(filename)

ax = axes[0]
M = results['M']
tlist = np.array(results['t_list'])
clustering_1 = results['clustering_1']
clustering_2 = results['clustering_2']
ax.plot(tlist/M, clustering_1,color=color_list[0],label='Original')
ax.plot(tlist/M, clustering_2,color=color_list[1],label='Stub matching')
ax.vlines(1, ymin=0,ymax=max(clustering_1), ls='--', colors=['grey'])
ax.legend(frameon=False)
ax.set_ylabel(r'$C_\mathrm{b}$')
ax.set_xlabel(r'Number of steps per edge')
ax.text(0., 1.1, f'(a) Lyon primary school', ha='left',va='center',
        transform=ax.transAxes,fontsize=font_size)


with open('../simulation/dat/coauth_mixing.pk','rb') as filename:
    results = pickle.load(filename)

ax = axes[1]
M = results['M']
tlist = np.array(results['t_list'])
clustering_1 = results['clustering_1']
clustering_2 = results['clustering_2']
ax.plot(tlist/M, clustering_1,color=color_list[0],label='Original')
ax.plot(tlist/M, clustering_2,color=color_list[1],label='Stub matching')
ax.vlines(1, ymin=0,ymax=max(clustering_1), ls='--', colors=['grey'])
ax.legend(frameon=False)
ax.set_ylabel(r'$C_\mathrm{b}$')
ax.set_xlabel(r'Number of steps per edge')
ax.text(0., 1.1, f'(b) Coauthorship in computer science', ha='left',va='center',
        transform=ax.transAxes,fontsize=font_size)



plt.show()
