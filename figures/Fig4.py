from gcm import *
from scipy.special import gamma
import matplotlib.pyplot as plt

#plot parameters
font_size=8
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font',family='serif',serif='Computer Modern')
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)
plt.rc('axes', labelsize=font_size)
colorlist = ["#c7e9b4","#41b6c4","#225ea8"]
width = 7.057/2
height = 1.2*width

fig, axes = plt.subplots(2,1, figsize=(width, height))
plt.subplots_adjust(left=0.3, bottom=0.1, right=0.87, top=0.92, wspace=0.18,
                    hspace=0.5)


#structure
mmin = 0
nmax = 10
mmax = 10
pn = np.zeros(nmax+1)
pn[4] = 0.5
pn[5] = 0.5
pn /= np.sum(pn)
gm = np.zeros(nmax+1)
gm[2] = 0.5
gm[3] = 0.5
gm /= np.sum(gm)

state_meta = get_state_meta(mmax, nmax, gm, pn)

#infection function
beta = lambda n,i,trate,nu: trate*i**nu
nu_c = bistability_threshold(beta,gm,pn, initial_params=(0.1,1.9))
print(nu_c)

#get map around invasion threshold
nu = 1.
trate_c = invasion_threshold(beta,gm,pn,fixed_args=(nu,))
print(trate_c)
trate_list = [0.05,trate_c,0.2]
label_list = [r"$\lambda < \lambda_c$", r"$\lambda = \lambda_c$", r"$\lambda > \lambda_c$"]

r_list = np.linspace(0.,0.3,100)
for trate,color,label in zip(reversed(trate_list),reversed(colorlist),
                          reversed(label_list)):
    inf_mat = infection_matrix(beta, nmax, args=(trate,nu))
    mf = lambda r: mf_map(r, inf_mat, state_meta)
    map_result = np.array([mf(r) for r in r_list])
    axes[0].plot(r_list, map_result, color=color,
                label=label, lw=1.5)

axes[0].plot(r_list, r_list, '--', color='black', lw=1.5)
axes[0].legend(frameon=False)
axes[0].set_xlabel(r"$r$")
axes[0].set_ylabel(r"$\mathcal{M}[\rho(r)]$")
axes[0].text(0.,1.1, '(a) Critical points', transform=axes[0].transAxes,
            fontsize = font_size+2)



#get map around trate_c for different mu
nu_list = [1.,nu_c,2.5]
label_list = [r"$\nu < \nu_c$", r"$\nu = \nu_c$", r"$\nu > \nu_c$"]
r_list = np.linspace(0.,0.8,100)
for nu,color,label in zip(reversed(nu_list),reversed(colorlist),
                          reversed(label_list)):
    trate = invasion_threshold(beta,gm,pn,fixed_args=(nu,))
    inf_mat = infection_matrix(beta, nmax, args=(trate,nu))
    mf = lambda r: mf_map(r, inf_mat, state_meta)
    map_result = np.array([mf(r) for r in r_list])
    axes[1].plot(r_list, map_result, color=color,
                label=label, lw=1.5)

axes[1].plot(r_list, r_list, '--', color='black', lw=1.5)
axes[1].legend(frameon=False)
axes[1].set_xlabel(r"$r$")
axes[1].set_ylabel(r"$\mathcal{M}[\rho(r)]$")
axes[1].text(0.,1.1, '(b) Tricritical points', transform=axes[1].transAxes,
            fontsize = font_size+2)
plt.savefig('../../paper/figs/Fig4.pdf')
# plt.show()

