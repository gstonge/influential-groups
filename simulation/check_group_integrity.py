import numpy as np
import pickle

#load dataset
dataset_dir = 'socio_data/'
dataset = 'coauth-dblp_sub'

#load graph data
with open(f'{dataset_dir}form_{dataset}.pk','rb') as filename:
    data = pickle.load(filename)

group_list = data['group_list']

max_label = 0
for i,group in enumerate(group_list):
    #check if same label
    set_group = set(group)
    if len(set_group) != len(group):
        print(f"multiplicity in group {i}")
        print(group)
