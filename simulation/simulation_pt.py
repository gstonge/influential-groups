import pickle
import json
import argparse
import numpy as np
from gcm import *
from _schon import PowerlawGroupSIS
import matplotlib.pyplot as plt
import horgg

#infection function
beta = lambda n,i,scale,shape: scale*i**shape

def get_bounds(beta,nmax,scale,shape):
    max_inf = nmax*beta(nmax,nmax,scale,shape) #upper bound
    min_inf = beta(1,1,scale,shape) #lower bound
    upper = max((max_inf,1))
    lower = min((min_inf,1))
    return (lower,upper)

def experience(conf, initial_infected_fraction, quasi=True):
    #load dataset
    dataset_dir = conf['dataset_dir'] #'socio_data/'
    dataset = conf['dataset']

    #load graph data
    with open(f'{dataset_dir}form_{dataset}.pk','rb') as filename:
        data = pickle.load(filename)

    seed = conf['seed']
    if 'nb_graph' in conf:
        nb_graph = conf['nb_graph']
    else:
        nb_graph = None

    gm = data['gm']
    pn = data['pn']
    nmax = data['nmax']
    orig_edge_list = data['edge_list']

    if nb_graph is None:
        edge_list_list = [orig_edge_list]
    else:
        if 'expansion_factor' in conf:
            factor = int(conf['expansion_factor'])
        else:
            factor = 1
        m_list = list(data['m_list'])*factor
        n_list = list(data['n_list'])*factor
        print(f"original size : {len(data['m_list'])}, new size: {len(m_list)}")
        #get new edge lists
        graph_generator = horgg.BCMS(m_list,n_list)
        horgg.BCMS.seed(seed) #optional, if nothing is given, it is seeded with time

        edge_list_list = [graph_generator.get_graph(mcmc_step=len(orig_edge_list))
                              for i in range(nb_graph)]

    #get invasion threshold and bistability threshold
    shape = conf['shape_infection']
    scale_c = invasion_threshold_safe(beta, gm, pn, fixed_args=(shape,),
                                      min_param=10**(-14), max_param=1)
    shape_c = bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                                         max_params=(1,7))
    print(f"Invasion threshold {scale_c}")
    print(f"Bistability threshold {shape_c}")

    #define simulation parameters and output
    lower_param = conf['lower_param']*scale_c
    higher_param = conf['higher_param']*scale_c
    nb_points = conf['nb_points']
    recovery_rate = conf['recovery_rate']
    nb_history = conf['nb_history']

    param_list = np.linspace(lower_param,higher_param,nb_points)
    mean_I_list = []
    std_I_list = []

    for scale in param_list:
        print(f"---------------------")
        print(f"scale parameter : {scale}")
        #varying parameters
        rate_bounds = get_bounds(beta,nmax,scale,shape)
        I_list = []
        for i,edge_list in enumerate(edge_list_list):
            print(f"simulation for network {i+1}...")
            #define process
            cont = PowerlawGroupSIS(edge_list,recovery_rate,scale,
                                    shape,rate_bounds)
            cont.infect_fraction(initial_infected_fraction)
            cont.seed(seed+i) #optional
            cont.initialize_history(nb_history)

            #define some measures
            cont.measure_prevalence()

            #evolve in the quasistationary state without measuring (burn-in)
            dt = conf['burnin_dt']
            dec_dt = conf['burnin_dec_dt']
            cont.evolve(dt,dec_dt,measure=False,quasistationary=quasi)

            #evolve and measure
            dt = conf['measure_dt']
            dec_dt = conf['measure_dec_dt']
            cont.evolve(dt,dec_dt,measure=True,quasistationary=quasi)

            #get the measure
            for measure in cont.get_measure_vector():
                name = measure.get_name()
                if name == "prevalence":
                    I_list.extend(measure.get_result())

        mean_I_list.append(np.mean(I_list))
        std_I_list.append(np.std(I_list))
        print(f"mean prevalence : {np.mean(I_list)}")
    return param_list, mean_I_list, std_I_list

#Define the parser
parser = argparse.ArgumentParser()
parser.add_argument("--dir", type=str, help="Experience directory")
parser.add_argument("--conf", type=str, help="Configuration file name")
args = parser.parse_args()

#load configuration file
with open(f"{args.dir}/{args.conf}.json") as filename:
    conf = json.load(filename)

#get the upper and lower branches
results = dict()

if conf['upper_branch']:
    print("========================")
    print("UPPER BRANCH")
    print("========================")
    param_list, mean_I_list_upper,std_I_list_upper = experience(
        conf, conf['init_I_upper'],conf['quasi_upper'])
    plt.plot(param_list,mean_I_list_upper)
    results['param_list'] = param_list
    results['mean_I_list_upper'] = mean_I_list_upper
    results['std_I_list_upper'] = std_I_list_upper

if conf['lower_branch']:
    print("")
    print("========================")
    print("LOWER BRANCH")
    print("========================")
    param_list, mean_I_list_lower,std_I_list_lower = experience(
        conf, conf['init_I_lower'], conf['quasi_lower'])
    plt.plot(param_list,mean_I_list_lower, '--')
    results['param_list'] = param_list
    results['mean_I_list_lower'] = mean_I_list_lower
    results['std_I_list_lower'] = std_I_list_lower

#plot
plt.show()

with open(f"{args.dir}/{args.conf}.pk", 'wb') as filename:
    pickle.dump(results,filename)

