import pickle
import networkx as nx
from networkx.algorithms import bipartite
from netwulf import visualize


# with open(f'./socio_data/form_thr3_LyonSchool.pk','rb') as filename:
with open(f'./socio_data/form_coauth-dblp_sub.pk','rb') as filename:
    data = pickle.load(filename)

edge_list = data['edge_list']

nodes = {f"n{u}" for u,v in edge_list}
groups = {f"g{v}" for u,v in edge_list}
edge_list = [(f"n{u}",f"g{v}") for u,v in edge_list]
print(len(nodes))

B = nx.Graph()
B.add_nodes_from(nodes, bipartite=0)
B.add_nodes_from(groups, bipartite=1)
B.add_edges_from(edge_list)
print(nx.is_connected(B))

visualize(B)
