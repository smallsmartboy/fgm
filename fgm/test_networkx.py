import networkx as nx
import community
import matplotlib.pyplot as plt

# G = nx.random_graphs.powerlaw_cluster_graph(3000, 1, .4)
#
path = "E:\\data\\graph\\lan_noattack\\20151102_gt_layer_6.csv"
g = nx.read_edgelist(path, delimiter=',', nodetype=str)
part = community.best_partition(g)

part = community.best_partition(g)
values = [part.get(node) for node in g.nodes()]

nx.draw_spring(g, cmap=plt.get_cmap('jet'), node_color=values, node_size=30, with_labels=False)

mod = community.modularity(part, g)
print("modularity:", mod)
plt.show()

