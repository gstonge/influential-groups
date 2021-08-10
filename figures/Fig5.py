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
colorlist = ["#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8", "#0c2c84"]
width = 7.057
height = width/4

fig, axes = plt.subplots(1,3, figsize=(width, height))
plt.subplots_adjust(left=0.07, bottom=0.23, right=0.98, top=0.97, wspace=0.30,
                    hspace=0.25)


#Left side figure
#=====================
with open('./dat/Fig3_bif2.pk', 'rb') as datafile:
    results = pickle.load(datafile)

for key in reversed(["struct1","struct2","struct3"]):
    structure = results[key]
    if key == "struct1":
        # label = r"$g_m = \delta_{m,6},\; p_n = \delta_{n,3}$"
        label = r"$p_n = \delta_{n,3}$"
        color = colorlist[0]
    elif key == "struct2":
        # label = r"$g_m = \delta_{m,4},\; p_n = \delta_{n,4}$"
        label = r"$p_n = \delta_{n,4}$"
        color = colorlist[2]
    elif key == "struct3":
        # label = r"$g_m = \delta_{m,3},\; p_n = \delta_{n,5}$"
        label = r"$p_n = \delta_{n,5}$"
        color = colorlist[4]

    for stability,solution in structure.items():
        if stability == "stable":
            axes[0].plot(solution['param_list'],
                         solution['infected_fraction_list'],'-', lw=1.5,
                         color=color,label=label)
        elif stability == "unstable":
            axes[0].plot(solution['param_list'],
                         solution['infected_fraction_list'],'--', lw=1.5,
                         color=color)
            axes[0].plot((0.96,solution['param_list'][-1]),
                         (0,0),'-', lw=1.5,
                         color=color)
axes[0].set_xticks([0.95,1.,1.05,1.1])
axes[0].set_xlabel(r"$\lambda/\lambda_\mathrm{c}$")
axes[0].set_ylabel(r"Infected fraction $I^*$")
xvals = [1.,1.03,1.06]
labelLines(list(axes[0].get_lines()), zorder=2.5,align=True,xvals=xvals,
           color='black',fontsize=font_size)
axes[0].text(0.12, 0.85, r'(a)', ha='center',
             va='center', transform=axes[0].transAxes,
            fontsize=font_size)


#middle side figure
#=====================
#infection params
nu_list = [1.5, 1.7, 1.9080447119702004, 2.1, 2.3]
nu_c = 1.9080447119702004
param_list = np.linspace(0.8,1.1,2000) #normalized with threshold
beta = lambda n,i,trate,nu: trate*i**nu

with open('./dat/Fig3_bif1.pk', 'rb') as datafile:
    results = pickle.load(datafile)

label_below_exist = False
label_above_exist = False
label_equal_exist = False
for j,nu in enumerate(reversed(nu_list)):
    color = colorlist[4-j]
    label = None
    if j in [1,2,3]:
        if nu < nu_c:
            label=r"$\nu < \nu_\mathrm{c}$"
        elif nu == nu_c:
            label=r"$\nu = \nu_\mathrm{c}$"
        elif nu > nu_c:
            label=r"$\nu > \nu_\mathrm{c}$"

    if nu == 2.1:
        axes[1].plot(results[nu]['stable']['param_list'][1:],
                 results[nu]['stable']['infected_fraction_list'][1:],
                 '-', lw=1.5, color=color,
                 label=label)
    else:
        axes[1].plot(results[nu]['stable']['param_list'],
                 results[nu]['stable']['infected_fraction_list'], '-',
                     lw=1.5, color=color, label=label,
                     zorder=(2 if nu != nu_c else -2))

    if nu > nu_c:
        axes[1].plot(results[nu]['unstable']['param_list'],
                 results[nu]['unstable']['infected_fraction_list'],'--',
                 lw=1.5, color=color)
        axes[1].plot((0.940,results[nu]['unstable']['param_list'][-1]),
                 (0,0),'-',
                 lw=1.5, color="#225ea8")
axes[1].set_xticks([0.95,1.,1.05,1.1])
axes[1].set_yticks([0.,0.2,0.4,0.6])
axes[1].set_xlabel(r"$\lambda/\lambda_\mathrm{c}$")
axes[1].set_ylabel(r"Infected fraction $I^*$")
xvals = [1.,1.03,1.06]
labelLines(list(axes[1].get_lines()), zorder=2.5,align=True,xvals=xvals,
           color='black',fontsize=font_size)
axes[1].text(0.12, 0.85, r'(b)', ha='center',
             va='center', transform=axes[1].transAxes,
            fontsize=font_size)

#right panel -- phase diagram
#--------------------------

#invasion threshold
#=====================
with open('./dat/Fig3_pd_inv.pk', 'rb') as datafile:
    data_inv = pickle.load(datafile)

axes[2].plot(np.array(data_inv['lambda_list'])/np.array(data_inv['lambda_list']),
         data_inv['nu_list'], '-', color='#1a1a1a', lw=2.)

#tricritical point
#===================
#membership
mmax = 3
gm = np.zeros(mmax+1)
gm[mmax] += 1

#group distribution
nmax = 4
pn = np.zeros(nmax+1)
pn[nmax] += 1

#infection
beta = lambda n,i,trate,nu: trate*i**nu
nu_c = bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                              max_params=(1,7))
lambda_c = invasion_threshold_safe(beta, gm, pn, fixed_args=(nu_c,),
                            min_param=10**(-14), max_param=1)
tricrit = (lambda_c,nu_c)
axes[2].scatter([tricrit[0]/lambda_c],[tricrit[1]], 80, marker='*',
                color='#1a1a1a',zorder=3)


#persistence threshold
#=====================
with open('./dat/Fig3_pd_pers.pk', 'rb') as datafile:
    data_pers = pickle.load(datafile)

#find index above tricritical point
rescaled_lambda = np.array(data_pers['lambda_list'])/np.array(
    data_inv['lambda_list'])
axes[2].plot(rescaled_lambda, data_pers['nu_list'], '--', color='#1a1a1a')
axes[2].text(0.12, 0.85, r'(c)', ha='center',
             va='center', transform=axes[2].transAxes,
            fontsize=font_size)
axes[2].set_xlabel(r"$\lambda/\lambda_\mathrm{c}$")
axes[2].set_ylabel(r"Nonlinear exponent $\nu$")
axes[2].set_xlim([0.,1.5])
axes[2].set_ylim([0.1,4])

#healthy region
axes[2].fill_between([0.,1.01*min(rescaled_lambda)], 4., color="#c7e9b4")
axes[2].fill_between(rescaled_lambda, data_pers['nu_list'] , color="#c7e9b4")
axes[2].text(0.32, 0.4, r'Healthy', ha='center',
             va='center', transform=axes[2].transAxes,
            fontsize=font_size)

#bistable region
axes[2].fill_between(rescaled_lambda, 4, data_pers['nu_list'] , color="#7fcdbb")
axes[2].text(0.52, 0.90, r'Bistable', ha='center',
             va='center', transform=axes[2].transAxes,
            fontsize=font_size)

#endemic
axes[2].fill_between([1,1.5], 4, color="#1d91c0")
axes[2].text(0.84, 0.40, r'Endemic', ha='center',
             va='center', transform=axes[2].transAxes,
            fontsize=font_size)

plt.show()


