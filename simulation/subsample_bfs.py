import pickle
import numpy as np
from collections import deque
from format_data import format_data

def bfs(adj_dict, group_list, source_group, depth_limit):
    observed_group = {source_group}
    depth = 0
    queue = deque([(source_group,depth)])
    while queue and queue[0][1] < depth_limit:
        group_parent, depth = queue.popleft()
        print(f"Group {group_parent}, depth {depth}")
        for node in group_list[group_parent]:
            for group_child in adj_dict[node]:
                if group_child not in observed_group:
                    observed_group.add(group_child)
                    queue.append((group_child,depth +1))
    return observed_group

datafile = './socio_data/form_coauth-dblp.pk'
new_datafile = './socio_data/form_coauth-dblp_sub.pk'

seed = 42
depth_limit = 3
np.random.seed(seed)

with open(datafile,'rb') as filename:
    data = pickle.load(filename)

group_list = data['group_list']
adj_dict = data['adj_dict']

#bfs subsample
source_group = np.random.randint(0,len(group_list))
observed_group = bfs(adj_dict, group_list, source_group, depth_limit)
print("Number of groups : ",len(observed_group))
subgroup_list = [group_list[group] for group in observed_group]
results = format_data(subgroup_list)

with open(new_datafile, 'wb') as output:
    pickle.dump(results,output)
