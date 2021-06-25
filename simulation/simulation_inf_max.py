import pickle
import json
import argparse
import numpy as np
from gcm import *
from _schon import PowerlawGroupSIS
import matplotlib.pyplot as plt
import horgg
import heapq

#infection function
beta = lambda n,i,scale,shape: scale*i**shape

def get_bounds(beta,nmax,scale,shape):
    max_inf = nmax*beta(nmax,nmax,scale,shape) #upper bound
    min_inf = beta(1,1,scale,shape) #lower bound
    upper = max((max_inf,1))
    lower = min((min_inf,1))
    return (lower,upper)

def influential_spreaders_strategy(inf_mat, data, nb_node):
    adj_dict = data['adj_dict']
    #get the maximal degree nodes in order
    node_list = np.array([node for node in adj_dict])
    m_list = np.array([len(adj_dict[node]) for node in adj_dict])
    indices = np.argsort(m_list)
    #sort the nodes
    node_list = list(reversed(node_list[indices]))
    Inode_set = set()
    for i in range(nb_node):
        Inode_set.add(node_list[i])
    return Inode_set

def random_strategy(inf_mat, data, nb_node):
    adj_dict = data['adj_dict']
    #get the maximal degree nodes in order
    node_list = np.array([node for node in adj_dict])
    Inode_set = set()
    Inode_set = set(np.random.choice(node_list,size=nb_node,replace=False))
    return Inode_set

def influential_groups_strategy(inf_mat, data, nb_node):
    #for each group assign cost-effective ratios Rni for each configuration
    #assumes a small number of nodes to infect---biased for high density
    group_list = data['group_list']
    Q = []
    heapq.heapify(Q)
    counter = 0
    for g, group in enumerate(group_list):
        n = len(group)
        for i in range(1,n-1):
            Rni = inf_mat[n,i]*(n-i)/i
            heapq.heappush(Q, (-Rni, counter, (n,i), g))
            counter += 1
    Inode_set = set()
    while len(Inode_set) < nb_node:
        item = heapq.heappop(Q)
        n,i = item[2]
        g = item[3]
        if (nb_node - len(Inode_set)) >= i:
            new_nodes = np.random.choice(group_list[g], size=i, replace=False)
        else:
            new_nodes = np.random.choice(group_list[g],
                                         size=nb_node-len(Inode_set),
                                         replace=False)
        for node in new_nodes:
            Inode_set.add(node)
    return Inode_set

strategy_map = {'random':random_strategy,
                'spreaders':influential_spreaders_strategy,
                'groups':influential_groups_strategy}

def experience(conf):
    #load dataset
    dataset_dir = conf['dataset_dir'] #'socio_data/'
    dataset = conf['dataset']

    #load graph data
    with open(f'{dataset_dir}form_{dataset}.pk','rb') as filename:
        data = pickle.load(filename)

    seed = conf['seed']
    gm = data['gm']
    pn = data['pn']
    nmax = data['nmax']
    edge_list = data['edge_list']

    #get invasion threshold and bistability threshold
    shape = conf['shape_infection']
    scale_c = invasion_threshold_safe(beta, gm, pn, fixed_args=(shape,),
                                      min_param=10**(-14), max_param=1)
    shape_c = bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14),1),
                                         max_params=(1,7))
    print(f"Invasion threshold {scale_c}")
    print(f"Bistability threshold {shape_c}")

    #define simulation parameters and output
    scale = conf['param']*scale_c
    recovery_rate = conf['recovery_rate']
    rate_bounds = get_bounds(beta,nmax,scale,shape)
    inf_mat = infection_matrix(beta, nmax, args=(scale,shape))

    #define contagion process
    I_vec_list = []
    # cont = PowerlawGroupSIS(edge_list,recovery_rate,scale,
                            # shape,rate_bounds)
    # cont.seed(seed) #optional
    # cont.measure_prevalence()

    #pick influence maximization strategy
    strategy = strategy_map[conf['strategy']]
    nb_node = conf['nb_node']
    Inode_set = strategy(inf_mat, data, nb_node)
    print(len(Inode_set))

    for j in range(conf['nb_runs']):
        print(f"simulation run {j+1}...")
        cont = PowerlawGroupSIS(edge_list,recovery_rate,scale,
                                shape,rate_bounds)
        cont.seed(seed+j) #optional
        cont.measure_prevalence()
        #reset and infect nodes
        cont.infect_node_set(Inode_set)

        #evolve and measure
        dt = conf['measure_dt']
        dec_dt = conf['measure_dec_dt']
        cont.evolve(dt,dec_dt,measure=True,quasistationary=False)

        #get the measure
        for measure in cont.get_measure_vector():
            name = measure.get_name()
            if name == "prevalence":
                I_vec_list.append(measure.get_result())

    return I_vec_list

#Define the parser
parser = argparse.ArgumentParser()
parser.add_argument("--dir", type=str, help="Experience directory")
parser.add_argument("--conf", type=str, help="Configuration file name")
args = parser.parse_args()

#load configuration file
with open(f"{args.dir}/{args.conf}.json") as filename:
    conf = json.load(filename)

#get the upper and lower branches
I_vec_list = experience(conf)
for Ivec in I_vec_list:
    print(len(Ivec))
I_vec_list = np.array(I_vec_list)
dt = conf['measure_dec_dt']
tvec = [0]
for i in range(len(I_vec_list[0])-1):
    tvec.append(tvec[-1]+dt)
assert len(tvec) == len(I_vec_list[0])

results = dict()
results['I_vec_list'] = I_vec_list
results['I_vec_mean'] = np.mean(I_vec_list, axis=0)
results['I_vec_std'] = np.std(I_vec_list, axis=0)
results['tvec'] = tvec

for Ivec in I_vec_list:
    plt.plot(tvec,Ivec, color='black', alpha=0.7)
plt.errorbar(tvec, results['I_vec_mean'], yerr=results['I_vec_std'],
             fmt='-', color='red')

#plot
plt.show()

with open(f"{args.dir}/{args.conf}.pk", 'wb') as filename:
    pickle.dump(results,filename)

