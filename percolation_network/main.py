import numpy as np
import matplotlib.pyplot as plt


class PercolationNetwork:
    def __init__(self, L):
        """
        Define a network object with lattice of length L x L
        """
        self.L = L
        self.lattice = [[[] for i in range(self.L)] for j in range(self.L)]
        self.AddedNodes = []
        self.Clusters = []


    def GetNeighbors(self, node):
        """
        Return the neighboring nodes of a given node
        """

        if node[0] == 0:
            left_bound = self.L - 1
            right_bound = 1
        elif node[0] == (self.L - 1):
            right_bound = -(self.L - 1)
            left_bound = -1
        else:
            left_bound = -1
            right_bound = 1
        
        if node[1] == 0:
            lower_bound = self.L - 1
            upper_bound = 1
        elif node[1] == (self.L - 1):
            upper_bound = -(self.L - 1)
            lower_bound = -1
        else:
            upper_bound = 1
            lower_bound = -1

        n1 = [node[0], node[1] + upper_bound]
        n2 = [node[0] + right_bound, node[1]]
        n3 = [node[0], node[1] + lower_bound]
        n4 = [node[0] + left_bound, node[1]]

        return [n1, n2, n3, n4]


    def AddEdge(self, node, p):
        neighbors = self.GetNeighbors(node)

        for i in neighbors:
            num = i[0] * self.L + i[1]
            if num in self.AddedNodes:
                neighbors.remove(i)
        for n in neighbors:
            if np.random.random() < p:
                self.lattice[node[0]][node[1]].append(n)
                self.lattice[n[0]][n[1]].append(node)
        self.AddedNodes.append(node[0]*self.L + node[1])


    def Display(self):
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_title("2D regular Percolation Network")
        for i in range(self.L):
            for j in range(self.L):
                ax.plot(i, j, marker='o', color='b')
                for n in self.lattice[i][j]:
                    if abs(n[0] - i) > 1 and not abs(n[1] - j) > 1:
                        ax.plot([i, i - 1/2], [j, n[1]], 'r-')
                        ax.plot([n[0], n[0] + 1/2], [j, n[1]], 'r-')
                    elif abs(n[1] - j) > 1 and not abs(n[0] - i) > 1:
                        ax.plot([i, n[0]], [j, j - 1/2], 'r-')
                        ax.plot([i, n[0]], [n[1], n[1] + 1/2], 'r-')
                    elif abs(n[1] - j) > 1 and abs(n[0] - i) > 1:
                        ax.plot([i, i - 1/2], [j, j - 1/2], 'r-')
                        ax.plot([n[0], n[0] + 1/2], [n[1], n[1] + 1/2], 'r-')
                    else:
                        ax.plot([i, n[0]], [j, n[1]], 'r-')

        plt.show()


    def FindClusterFromNode(self, node, visited):
        cluster = []
        if (node[0] * self.L + node[1]) not in visited:
            cluster.append(node)
            current_neighbor = self.lattice[node[0]][node[1]]
            while len(current_neighbor) != 0:
                next_neighbor = []
                for n in current_neighbor:
                    if (n[0] * self.L + n[1]) not in visited:
                        cluster.append(n)
                        visited.append(n[0] * self.L + n[1])
                        for node in self.lattice[n[0]][n[1]]:
                            if (node[0] * self.L + node[1]) not in visited:
                                next_neighbor.append(node)
                current_neighbor = next_neighbor
        return cluster, visited




    def FindAllClusters(self):
        visited = []
        for i in range(self.L):
            for j in range(self.L):
                cluster, visited = self.FindClusterFromNode([i, j], visited)
                if len(cluster) > 0:
                    self.Clusters.append(cluster)
        



L = 10
regular = PercolationNetwork(L)

for i in range(L):
    for j in range(L):
        regular.AddEdge([i, j], 0.4)

regular.FindAllClusters()
print(len(regular.Clusters))

regular.Display()