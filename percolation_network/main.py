import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, ax = plt.subplots(2, 2, figsize=(10, 10))


def grid_percolation(m, n, p, periodic=True, create_using=None):
    """Returns the two-dimensional percolated grid graph
    
    The percolated grid graph disconnected edges by probability p

    Parameters
    ----------
    m, n : int or iterable container of nodes
        If an integer, nodes are from `range(n)`.
        If a container, elements become the coordinate of the nodes.
    
    p : int
        Each connected edges from the lattice has probability p to be removed.

    periodic : bool or iterable
        If `periodic` is True, both dimensions are periodic. If False, none
        are periodic.  If `periodic` is iterable, it should yield 2 bool
        values indicating whether the 1st and 2nd axes, respectively, are
        periodic.

    create_using : NetworkX graph constructor, optional (default=nx.Graph)
        Graph type to create. If graph instance, then cleared before populated.
    """
    G = nx.grid_2d_graph(m, n, periodic, create_using)
    for e in G.edges:
        if np.random.random() < p:
            n1, n2 = e[0], e[1]
            G.remove_edge(n1, n2)
    return G


def plot_grid_clusters(m, n, p):
    G = grid_percolation(m, n, p)
    clusters = list(nx.connected_components(G))
    indices = np.argsort([len(x) for x in clusters])
    sorted_indices = indices[::-1]
    sorted_clusters = np.array(clusters)[sorted_indices]
    color_list = []

    for i in range(len(sorted_clusters)):
        r, g, b = (i/len(clusters))**3, np.abs(np.sin(np.pi*i/2/len(clusters))), np.random.random()
        color_list.append((r, g, b))

    for i, cluster in enumerate(sorted_clusters):
        color = color_list[i]
        visted_nodes = []
        for c in cluster:
            visted_nodes.append(c)
            ax[0][0].plot(c[0], c[1], "o", color=color, markersize=4)
            for e in G.neighbors(c):
                if e not in visted_nodes:
                    x = []
                    y = []
                    n1, n2 = c, e
                    # if top row connect with the bottom row
                    if (n1[1] == 0 and n2[1] == n-1) or (n2[1] == 0 and n1[1] == n-1):
                        x.append(n1[0])
                        y.append(n1[1])
                        x.append(n1[0])
                        y.append(n1[1] + 1/2)

                        x.append(n2[0])
                        y.append(n2[1])   
                        x.append(n2[0])
                        y.append(n2[1] - 1/2)
                    # if the left column connect with the right column
                    elif (n1[0] == 0 and n2[0] == m-1) or (n2[0] == 0 and n1[0] == m-1):
                        x.append(n1[0])
                        y.append(n1[1])
                        x.append(n1[0] + 1/2)
                        y.append(n1[1])

                        x.append(n2[0])
                        y.append(n2[1])
                        x.append(n2[0] - 1/2)
                        y.append(n2[1])
                    
                    else:
                        x.append(n1[0])
                        y.append(n1[1])
                        x.append(n2[0])
                        y.append(n2[1])
                
                if len(x) > 3:
                    ax[0][0].plot(x[:2], y[:2], "-", color=color)
                    ax[0][0].plot(x[2:], y[2:], "-", color=color)
                else:
                    ax[0][0].plot(x, y, "-", color=color)
        ax[0][0].axis("off")


def triangular_percolation(m, n, p, periodic=True, create_using=None):
    """Returns the triangular percolated graph
    
    The percolated graph disconnected edges by probability p

    Parameters
    ----------
    m, n : int or iterable container of nodes
        If an integer, nodes are from `range(n)`.
        If a container, elements become the coordinate of the nodes.
    
    p : int
        Each connected edges from the lattice has probability p
        to be removed.

    periodic : bool or iterable
        If `periodic` is True, both dimensions are periodic. If False, none
        are periodic.  If `periodic` is iterable, it should yield 2 bool
        values indicating whether the 1st and 2nd axes, respectively, are
        periodic.

    create_using : NetworkX graph constructor, optional (default=nx.Graph)
        Graph type to create. If graph instance, then cleared before populated.
    """
    G = nx.triangular_lattice_graph(m, n, periodic=False, create_using=None)
    for n in list(G.nodes):
        if np.random.random() < p:
            G.remove_node(n)
    return G


def plot_triangular_cluster(G, m, n):
    clusters = list(nx.connected_components(G))
    indices = np.argsort([len(x) for x in clusters])
    sorted_indices = indices[::-1]
    sorted_clusters = np.array(clusters)[sorted_indices]

    for i, c in enumerate(sorted_clusters):
        r, g, b = (i/len(clusters))**3, np.abs(np.sin(np.pi*i/2/len(clusters))), np.random.random()
        radius = 1/np.sqrt(3)
        for node in c:
            pos = nx.get_node_attributes(G, "pos")[node]
            hexagon = patches.RegularPolygon(pos, 6, radius=radius, facecolor=(r, g, b))
            ax[1][0].add_patch(hexagon)
    
    ax[1][0].set_aspect('equal')
    ax[1][0].set_xlim(-1, n/2+1.5)
    ax[1][0].set_ylim(-1, m-1/2)
    ax[1][0].axis("off")




m, n, p = 20, 20, 0.5
plot_grid_clusters(m, n, p)

m, n, p = 20, 40, 0.5

G = triangular_percolation(m, n, p)

plot_triangular_cluster(G, m, n)



m, n, p = 1024, 1024, 0.5

G = grid_percolation(m, n, p)

clusters = list(nx.connected_components(G))

max_val = 0
max_index = []
for i, j in enumerate(clusters):
    if len(j) > max_val:
        max_index = [i]
        max_val = len(j)
    elif len(j) == max_val:
        max_index.append(i)

print(max_val, max_index)

C = nx.Graph()
C.add_nodes_from(clusters[max_index[0]])

x, y = np.array(list(clusters[max_index[0]])).T[0], np.array(list(clusters[max_index[0]])).T[1]
ax[0][1].plot(x, y, 'o', markersize=1/np.log(m*n), color="k")
ax[0][1].axis("off")



m, n, p = 1024, 2*1024, 0.5

G = triangular_percolation(m, n, p)
clusters = list(nx.connected_components(G))


max_val = 0
max_index = []
for i, j in enumerate(clusters):
    if len(j) > max_val:
        max_index = [i]
        max_val = len(j)
    elif len(j) == max_val:
        max_index.append(i)

print(max_val, max_index)

C = nx.Graph()
C.add_nodes_from(clusters[max_index[0]])

x, y = np.array(list(clusters[max_index[0]])).T[0], np.array(list(clusters[max_index[0]])).T[1]
ax[1][1].plot(x, y, 'o', markersize=1/np.log(m*n), color="k")
ax[1][1].axis("off")

fig.tight_layout()

plt.savefig("./figures/universality.png")