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

fig, axes = plt.subplots(1,2, figsize=(width, height))
plt.subplots_adjust(left=0.17, bottom=0.22, right=0.88, top=0.98, wspace=0.45,
                    hspace=0.25)


#Left side figure
#=====================
with open('./dat/Fig5_bis_thresh.pk', 'rb') as datafile:
    results = pickle.load(datafile)

label_dict = {100:r"$10^2$", 300:r"$300$", 1000:r"$10^3$",
              3000:r"$3000$", 10000:r"$10^4$"}

first_element = next(results.__iter__())
gamma_list = list(results[first_element].keys())
mmax_list = [mmax for mmax in results]

for i,mmax in enumerate([100,1000,10000]):
    bthreshold = [results[mmax][gamma] for gamma in gamma_list]
    axes[1].plot(gamma_list,bthreshold,'-', lw=1.5,
                     color=colorlist[2*i],
                 label=r"$m_\mathrm{max}$="+label_dict[mmax])#+f"{mmax}")

leg = axes[1].legend()
leg.get_frame().set_linewidth(0.0)
axes[1].set_xlabel(r"Membership exponent $\gamma_m$")
axes[1].set_ylabel(r"Bistability threshold $\nu_\mathrm{c}$")

axes[1].text(0.12, 0.85, r'(b)', ha='center',
             va='center', transform=axes[1].transAxes,
             fontsize=font_size)

#Right side figure
#=====================
# plot_gamma = [2.1,3.,3.8,4.,4.2]
plot_gamma = [2.1,3.,4.2]
for i,gamma in enumerate(plot_gamma):
    bthreshold = [results[mmax][gamma] for mmax in mmax_list]
    axes[0].semilogx(mmax_list,bthreshold,'-', lw=1.5,
                     color=colorlist[2*i],
                 label=r"$\gamma_m$="+f"{gamma}")

leg = axes[0].legend(loc=(0.18,0.57))
leg.get_frame().set_linewidth(0.0)
axes[0].set_xlabel(r"Maximal membership $m_\mathrm{max}$")
axes[0].set_ylabel(r"Bistability threshold $\nu_\mathrm{c}$")

axes[0].text(0.12, 0.85, r'(a)', ha='center',
             va='center', transform=axes[0].transAxes,
             fontsize=font_size)



plt.savefig('../../paper/figs/Fig7.pdf')
plt.show()


