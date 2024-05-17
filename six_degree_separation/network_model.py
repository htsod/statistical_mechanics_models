import numpy as np
import matplotlib.pyplot as plt


class  network:
    def __init__(self, L, Z, p):
        """
        Network object with L nodes Z edges with a periodic boundary condition
        """
        self.nodes = [i for i in range(L)]
        self.L = L
        self.Z = Z
        self.p = p
        edges = [[] for i in range(L)]

        for i in range(L):
            for j in range(int(Z/2)):
                if (i + 1 + j) >= len(self.nodes):
                    upper_bound = i + 1 + j - len(self.nodes)
                    lower_bound = i - 1 - j
                else:
                    upper_bound = i + 1 + j
                    lower_bound = i - 1 - j
                edges[i].append(self.nodes[lower_bound])
                edges[i].append(self.nodes[upper_bound])
        self.edges= edges
    
    
    def HasNode(self, node):
        """
        Check if the network contains a node
        """
        if node in self.nodes:
            return True
        else:
            return False


    def AddNode(self, node):
        """
        Add a node to the network
        """
        if not (self.HasNode(node)):
            self.nodes.append(node)
            self.nodes.sort()
            node_index = self.nodes.index(node)
            self.edges.insert(node_index, [])
            for i in range(int(self.Z/2)):
                if (node_index + 1 + i) >= len(self.nodes):
                    lower_bound = node_index - 1 - i
                    upper_bound = node_index + 1 + i - len(self.nodes)
                else:
                    lower_bound = node_index - 1 - i
                    upper_bound = node_index + 1 + i

                self.edges[node_index].append(self.nodes[lower_bound])
                self.edges[node_index].append(self.nodes[upper_bound])

                self.edges[lower_bound].remove(self.nodes[upper_bound])
                self.edges[lower_bound].append(self.nodes[node_index])

                self.edges[upper_bound].remove(self.nodes[lower_bound])
                self.edges[upper_bound].append(self.nodes[node_index])
        else:
            print("The node already exist")

    
    def AddEdge(self, node1, node2):
        """
        Add an edge that connects node1 and node2
        """
        if self.HasNode(node1) and self.HasNode(node2):
            n1_index = self.nodes.index(node1)
            n2_index = self.nodes.index(node2)
            if (node2 in self.edges[n1_index]) or (node2 in self.edges[n2_index]):
                print("Given nodes already connected")
            else:
                self.edges[n1_index].append(node2)
                self.edges[n2_index].append(node1)
        else:
            print("One of the node does not exist")

    
    def RemoveEdge(self, node1, node2):
        if self.HasNode(node1) and self.HasNode(node2):
            n1_index = self.nodes.index(node1)
            n2_index = self.nodes.index(node2)
            if (node2 in self.edges[n1_index]) or (node2 in self.edges[n2_index]):
                self.edges[n1_index].remove(node2)
                self.edges[n2_index].remove(node1)
            else:
                print("The given nodes are not connected")
        else:
            print("One or more of the given nodes are not in the network")



    def RandomAdd(self):
        """
        Randomly add p*L*Z/2 random shortcuts to the network
        """
        iteration = int(self.L * self.p * self.Z/2)
        for i in range(iteration):
            n1 = np.random.choice(self.nodes)
            n2 = np.random.choice(self.nodes)
            if n1 != n2:
                self.AddEdge(n1, n2)
            else:
                print("The nodes are equal")


    def GetNodes(self):
        """
        Return the nodes of the network
        """
        return self.nodes

    
    def GetNeighbors(self, node):
        """
        Print the neighbors of the given node
        """
        if self.HasNode(node):
            node_index = self.nodes.index(node)
            return self.edges[node_index]
        else:
            print("The node does not exist")


    def _NodeCoordinate_(self):
        """
        Project the nodes and edges on a polar coordinate
        """
        radius = len(self.nodes) * 10
        ang_interval = 2 * np.pi / len(self.nodes)
        nodes_coor = {}
        for i in range(len(self.nodes)):
            x = radius * np.cos(i * ang_interval)
            y = radius * np.sin(i * ang_interval)
            nodes_coor.update({self.nodes[i]: [x, y]})
        edges = []
        for i in range(len(self.nodes)):
            x_data = []
            y_data = []
            for j in range(len(self.edges[i])):
                x_data.append(nodes_coor[self.edges[i][j]][0])
                y_data.append(nodes_coor[self.edges[i][j]][1])
            edges.append([x_data, y_data])
        return nodes_coor, edges
        

    def Display(self):
        """
        Display the network on a polar coordinate
        """
        nodes_coor, edges= self._NodeCoordinate_()
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_axis_off()
        for i in range(len(self.nodes)):
            ax.scatter(nodes_coor[self.nodes[i]][0], nodes_coor[self.nodes[i]][1], color="b", s=5*2**len(self.edges[i]))
            for j in range(len(edges[i][0])):
                x1 = [nodes_coor[self.nodes[i]][0], edges[i][0][j]]
                y1 = [nodes_coor[self.nodes[i]][1], edges[i][1][j]]
                ax.plot(x1, y1, color='red')
        plt.show()


    def FindPathLengthsFromNode(self, node):
        """
        Find the path length to every another node from the given node
        """
        currentShell = [node]
        distance = 1
        distance_list = [[node], [0]]
        while len(distance_list[0]) != len(self.nodes):
            nextShell = []
            for n in currentShell:
                for i in self.edges[n]:
                    if i not in distance_list[0]:
                        distance_list[0].append(i)
                        distance_list[1].append(distance)
                    if i not in nextShell:
                        nextShell.append(i)
            distance += 1
            currentShell = nextShell
        return distance_list[1]


    def FindAllPathLengths(self):
        """
        Find all path length of every nodes
        """
        distance_list = []
        for n in self.nodes:
            distances = self.FindPathLengthsFromNode(n)
            for i in distances:
                distance_list.append(i)
        return distance_list


    def WattsStrogaztzProcedure(self):
        """
        Construct a network using Watts Strogaztz Procedure
        """
        for node in self.nodes:
            node_index = self.nodes.index(node)

            for i in range(int(self.Z/2)):
                if (node_index + 1 + i) >= len(self.nodes):
                    lower_bound = node_index - 1 - i
                    upper_bound = node_index + 1 + i - len(self.nodes)
                else:
                    lower_bound = node_index - 1 - i
                    upper_bound = node_index + 1 + i
                if np.random.random() < self.p:
                    condition = 0
                    while condition < 1:
                        if (lower_bound in self.edges[node_index]) and (upper_bound in self.edges[node_index]):
                            self.edges[node_index].remove(lower_bound)
                            self.edges[node_index].remove(upper_bound)
                        node_add = np.random.choice(self.nodes)
                        if node != node_add and (node_add not in self.edges[node_index]):
                            self.AddEdge(node, node_add)
                            condition += 1


    def FindAveragePathLength(self):
        """
        Calculate the average of all path lemgths
        """
        distance_list = self.FindAllPathLengths()
        return distance_list, np.array(distance_list).mean()
    

    def FindBetweenessFromNode(self, node):
        """
        Calculate the betweeness from a single node
        """
        pass


    def FindAllBetweeness(self):
        """
        Find all betweeness of all nodes
        """
        pass
