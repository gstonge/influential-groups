import json
import pickle
import numpy as np
import itertools
from os import path

def format_data(group_list):
    #remove possible multiplicity in groups
    group_list_ = [list(set(g)) for g in group_list]

    #get group size list
    n_list = [len(g) for g in group_list_]

    #get node adjacency
    node_mapping = dict() #relabel nodes from 0 to N-1
    adj_dict = dict()
    for g_id,g in enumerate(group_list_):
        for node in g:
            if node in node_mapping:
                adj_dict[node_mapping[node]].append(g_id)
            else:
                node_id = len(node_mapping)
                node_mapping[node] = node_id
                adj_dict[node_id] = [g_id]

    #relabel group list
    new_group_list = [[node_mapping[node] for node in group]
                      for group in group_list_]

    #get membership
    m_list = [len(adj_dict[node]) for node in adj_dict]
    nmax = max(n_list)
    mmax = max(m_list)

    #get edge-list
    edge_list = []
    for i in adj_dict:
        for j in adj_dict[i]:
            edge_list.append((i,j))

    #get distributions
    pn = np.zeros(nmax+1)
    for n in n_list:
        pn[n] += 1
    pn /= np.sum(pn)
    gm = np.zeros(mmax+1)
    for m in m_list:
        gm[m] += 1
    gm /= np.sum(gm)

    #get joint and conditional distributions
    Pmn = np.zeros((nmax+1,mmax+1))
    for node in adj_dict:
        for g_id in adj_dict[node]:
            Pmn[n_list[g_id],m_list[node]] += 1
    Pmn /= np.sum(Pmn) #joint
    Pm_n = np.zeros((nmax+1,mmax+1))
    Pn_m = np.zeros((mmax+1,nmax+1))
    for n in range(nmax+1):
        norm = np.sum(Pm_n[n,:])
        if norm > 0:
            Pm_n[n,:] /= norm
    for m in range(mmax+1):
        norm = np.sum(Pn_m[m,:])
        if norm > 0:
            Pn_m[m,:] /= norm

    results = dict()
    results['group_list'] = new_group_list
    results['adj_dict'] = adj_dict
    results['edge_list'] = edge_list
    results['m_list'] = m_list
    results['n_list'] = n_list
    results['pn'] = pn
    results['gm'] = gm
    results['Pmn'] = Pmn
    results['Pm_n'] = Pm_n
    results['Pn_m'] = Pn_m
    results['nmax'] = nmax
    results['mmax'] = mmax

    return results

dataset_dir = 'socio_data/'
n_minutes = 15
dataset_list = ['InVS13','InVS15','LH10','LyonSchool','SFHH','Thiers13']
thr_list = [1,3,5]
for dataset,thr in itertools.product(dataset_list,thr_list):
    filename = dataset_dir+'aggr_'+str(n_minutes)+'min_cliques_thr'+str(thr)+'_'+dataset+'.json'
    if path.exists(filename):
        group_list = json.load(open(filename,'r'))
        results = format_data(group_list)
        new_filename = dataset_dir+'form_thr'+str(thr)+'_'+dataset+'.pk'
        with open(new_filename, 'wb') as output:
            pickle.dump(results,output)

#other datasets
dataset_list = ['email-Eu_simplices', 'coauth-DBLP_simplices']
outname_list = ['email-eu', 'coauth-dblp']

for dataset,outname in zip(dataset_list,outname_list):
    filename = dataset_dir+dataset+'.json'
    if path.exists(filename):
        group_list = json.load(open(filename,'r'))
        results = format_data(group_list)
        new_filename = dataset_dir+"form_"+outname+'.pk'
        with open(new_filename, 'wb') as output:
            pickle.dump(results,output)
