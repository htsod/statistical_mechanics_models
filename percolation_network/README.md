## Percolation Network

A __Percolation Network__ is a fascinating computational model that demonstrates the qualitative changes that occur as "holes" are introduced into a two-dimensional network. This exercise provides a rich perspective on percolation through the lens of __continuous phase transitions__.

Starting with a well-connected lattice network, where each node is connected to its nearest neighbors, the percolation process simulates the progressive removal of nodes or edges. The probability of a connection being removed is denoted as $p$. As $p$ crosses a critical threshold, $p_{c}$, the network undergoes a dramatic transformation from a single, well-connected structure to scattered, disjoint clusters. How can we understand this emergent behavior?

### Phase Transition

For $p < p_{c}$, the network remains connected, with a large cluster dominating the structure. In this regime, the system exhibits __broken symmetry__, with long-range order enabling traversal through the largest cluster.

As $p > p_{c}$, the network fragments into smaller clusters, resembling a gas-like state with __translation symmetry__. This marks the __symmetric phase__ of the system.

This transition from a symmetric to a broken-symmetric state as $p$ varies is an example of a __phase transition__ in __statistical mechanics__. It is analogous to the transition between ferromagnetic and paramagnetic states, as described by the Ising model. For a detailed exploration of phase transitions through the Ising model, visit[Solution to the Ising Model](https://htsod.github.io/posts/ising_model/).


### Edges Percolation

In edge percolation, connections (edges) are removed with probability $p$. Below is an example of edge percolation on a periodic $20\times 20$ lattice with $p=0.4$. 

![grid_percolation](/percolation_network/figures/grid_percolation.png)

Nodes of different colors represent different clusters. At $p=0.5$, the system appears ambiguous—neither fully connected nor entirely fragmented. However, a slight change in $p$ reveals stark differences in the network's state:


![grid_percolation_transition](/percolation_network/figures/grid_percolation_transition.png)

- __Left__: At $p=0.4$, below the transition point, a single giant cluster dominates.

- __Right__: At $p=0.6$, above the transition point, the network fragments into smaller clusters.

It could be shown that at the critical probability $p_{c} = 0.5$, the network transitions from a fully connected state to disjoint clusters. Suppose we construct a dual graph of the original one that has $1-p$ probability of breaking edges. Both graph will undergo transition at the same $p$ value because they are statistically identical. Then, $1-p = p$, which proves $p_{c}= 0.5$


### Bond Percolation

In contrast to edge percolation, bond percolation involves the removal of nodes. For example, on a triangular lattice, bond percolation can be represented as a hexagonal tiling:

![Bond Percolation](/percolation_network/figures/triangular_percolation.png)

As $p$ varies from $p = 0.4$ to $p=0.5$, the system again exhibits a phase transition. Both edge and bond percolation display the same qualitative behavior, exemplifying __universality__—the insensitivity of phase transition properties to microscopic details.


### Universality and Scale invariant

__Universality__ states that by knowing only the symmetries on either side of the phase transition and nothing whatsoever about the microscopic Hamiltonian, we can already deduce, real observable properties of the system near its phase transition. [[1]](#1)

This observation leads to a powerful classification of phase transition across different system with different microscopic detials. In the example of bond and edge percolation, they share the same symmetry changes, and thus will be classify into the same universality class.

The universality could also be visualized by comparing them in the large $N$ limits at the critical point $p_{c}= 0.5$. The largest cluster of each case have the same statistical properties:

![universality](/percolation_network/figures/universality.png)

Furthermore, these clusters demonstrate __scale invariance__ at the critical point, maintaining their statistical properties under rescaling, which is a common features of __continuous phase transition__ at the critical point.

To explore the concepts of universality and scale invariant on a firm quantiative ground, recall the relevent parameters and observables for the percolation network:

- $L$ the length of the square network 

- $p$ probability of removing a connection

- $S$ the cluster size of a network

- $n(S)$ distribution of the cluster size

- $P(p) = \frac{S_{largest}(p)}{L*L}$ the relative size of the largest cluster as a function of $p$

<!-- At the $p \rightarrow 0$ limit, the graphs are fully connected with the nearest neighbors, and the cluster size of the network will be the same as the number of nodes in the network. $S(p=0) = L\times L$ and hence $P(p)=1$ . On the opposite extreme, $p \rightarrow 1$, the graph falls into disconnected nodes with $P(p) = 1/(L\timesL)$. -->

At the critical point $p_{c}$, the distribution of the cluster will scale with the following relation $n(S) \propto S^{-\tau}$ where $\tau$ has the theoretical prediction of $187/91$. The exponent is known to be the scaling exponent and is the same in the same universality class. The constant of proportionality might differ, but the scaling exponents are the same for the bond percolation and site percolation. To see this:

![distribution of the cluster size](/percolation_network/figures/cluster_distribution.png)

The gradient of the log-log plot directly compares the exponents. The above graph demonstrate that though differs in microscopic details, the scaling relation holds.

Now consider for $p$ around the critical point, the theoretical prediction tells us that $P(p) \propto (p_{c} - p)^{\beta}$ for $\beta=5/36$. Plotting the result:

![scaling near pc](/percolation_network/figures/log_log_plot_scaling_relation_near_pc.png)

You might notice that in the above graph, the experimental results do not align well with the theoretical prediction. This is the scaling exponent only valids numerically at the thermodynamics limit. Experimentally, it is impossible to simulate at the thermodynamics limit at $\lim_{L\rightarow \infty}$, but we make use of the fact that as $L$ go larger, the scaling function will collapse to the thermodynamics limit. This is also known from the method of __Finite-size scaling__.

![finite-size scaling](/percolation_network/figures/collapse_into_ingle_scaling_function_bond.png)




### Reference

<a id="1">[1]</a> 
Beekman, Aron and Rademaker, Louk and van Wezel, Jasper (2019). 
An introduction to spontaneous symmetry breaking. 
Phase transitions, 5(5), 93-94.