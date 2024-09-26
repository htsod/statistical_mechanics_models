## Generating Random Walks

This section explores the statistical properties and emergent symmetries of random walks.

### Random walk distribution

In a random walk, each step is determined randomly within the constraints of the allowed dimensionality. The goal of the first experiment is to study the statistical properties—such as the mean and standard deviation—of the final positions of the walkers after taking $N$ steps. In the continuum limit of a discrete random walk, the behavior of the system transitions into a differential equation:

$$ \frac{\delta \rho}{\delta t} = D \nabla^{2} $$

This is recognized as the __diffusion equation__, where $D$ is the diffusion constant, and $\rho$ is the evolving density of the walkers. The solution to this equation, representing the distribution of walkers, follows a Gaussian distribution.

To test how well the Gaussian approximation applies to random walks with a small number of steps, we compare the distributions for step sizes of 1, 2, 3, and 5:


![one_d_random_walks](/generating_random_walks/figures/one_d_vary_steps.png)

We observe that with only one step, the distribution deviates significantly from the Gaussian. However, even with just two steps, the distribution begins to resemble a Gaussian curve. As the number of steps increases, the approximation improves.



### Emergent symmetry

If walkers take only one step, their positions are confined to a square lattice, restricting the range of possible locations. However, when additional steps are allowed, the distribution evolves into a circular pattern. This is due to the fluctuations in the walk, which allow walkers to explore neighboring regions, resulting in a configuration with __spherical symmetry__.

This emergent symmetry is demonstrated in the figure below:

![two_dimensional_scatter](/generating_random_walks/figures/two_dimensional_scatter.png)


The square at the center represents walkers after taking just one step, while the blue dots represent walkers that have taken ten steps. This phenomenon, where symmetry arises with an increasing number of steps or fluctuations, is an example of emergent symmetry.


### Possible extension

By analyzing the symmetry of a system's distribution, we can infer its underlying evolution mechanisms. In the case of random walks, if we assume the microscopic rules of movement are unknown and we only have snapshots of the system at different stages, we can still deduce key properties. Initially, the system forms a square distribution, which later transitions into a circular one as more steps are taken. This change suggests that the microscopic rules allow only two-fold states for each step.

If we limit the number of steps back to one, a __second-order phase transition__ occurs, where the circular symmetry breaks down into cubic symmetry. This is analogous to the two-state spin system in the Ising model, where the step size is replaced by the critical temperature.

In today's data-rich environment, we often only have access to the system's distribution in configuration space, while the underlying dynamics may be unknown. This problem illustrates that by studying the symmetry patterns over time, we can infer the microscopic behavior governing the system.