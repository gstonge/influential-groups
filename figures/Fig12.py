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
width = 7.057/2
height = width/2.
nb_dist = 6
fig, small_axes = plt.subplots(ncols=1, nrows=nb_dist, figsize=(width,height))

plt.subplots_adjust(left=0.25, bottom=0.17, right=0.85, top=1.15, wspace=0.37,
                    hspace=0.25)


with open('./dat/Fig11_temp_evo.pk', 'rb') as datafile:
    results = pickle.load(datafile)
fnivec = results['fni']

nmax = 20
tlist = [50,60,70,80,90,100]
fni_dict = {fr"$t = {t}$":{i/nmax: fnivec[t][nmax][i] for i in range(nmax+1)} for t in tlist}
cm = LinearSegmentedColormap.from_list('cmap', ["#7fcdbb","#0c2c84"])
color_map = [cm(x) for x in np.linspace(0.8,0.1,nb_dist)]

joyplot(small_axes,fni_dict, xmin=-0.2, xmax=1., smooth=False,rescale=True,
        xticks=[0,0.2,0.4,0.6,0.8], color_map=color_map,labelpos=(-0.22,0.02))

for ax in small_axes:
    ax.set_xticks([])
    l, b, w, h = ax.get_position().bounds
    ax.set_position([l, 0.09+0.65*b, w-0.01, h])
small_axes[-1].set_xlabel(r'Infected fraction $i/n$')
small_axes[-1].set_xticks([0,0.2,0.4,0.6,0.8,1])
small_axes[0].text(0.25,1.38,r'Distribution $f_{n_\mathrm{max},i}$',fontsize=font_size,
                      transform=small_axes[0].transAxes)

plt.show()


