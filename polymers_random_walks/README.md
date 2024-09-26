## Polymers Random Walks

Polymers consist of small building blocks called monomers, which are bonded together to form a one-dimensional chain. When placed in an aqueous environment, the configuration of a polymer at any given time resembles a __self-avoiding random walk__.

In this project, we explore a novel measure for random walks: the __root mean square (RMS)__ distance between the starting point and the end of the walk. Intuitively, as the number of steps (or monomers) increases, the RMS distance grows. Additionally, for an __unrestricted random walk__, the probability of retracing previous steps is higher, which results in a shorter RMS distance compared to self-avoiding walks. The central question we investigate is whether the two types of random walks—unrestricted and self-avoiding—are intrinsically different. Given the differences in microscopic evolution between these walks, do they scale similarly as the number of steps increases?


### Numerical simulation

To answer this question, we simulate self-avoiding random walks and compare their RMS distances to the theoretical RMS distance of an unrestricted random walk, which scales as $\sqrt{N}$, where $N$ is the number of steps.


![self_avoiding](/polymers_random_walks/figures/exponent.png)


The results demonstrate a clear qualitative difference between the two types of random walks. In technical terms, as the length scale increases, the fluctuations do not wash out the microscopic differences, placing them in different universality classes. Specifically, the self-avoiding random walk scales as $N^{\frac{3}{4}}$, slightly larger than the $N^{\frac{1}{2}}$ scaling of the unrestricted random walk.

### Implications

As we vary the length scale from the microscopic evolution to macroscopic observable, some of the asymetrical states in a given system will be filled, essentially becoming insignificant in the thermodynamics limit ($N\rightarrow \infty$). However, in self-avoiding random walks, the microscopic features persist across larger length scales, leading to significant differences in statistical properties.

From a theoretical standpoint, this project raises the question of which microscopic details persist and which are "erased" as the system grows. If we figure out this underlying truth, we could directly assign the universality class from a system microscopic degree of freedom and predicts its statistical macroscopic observables directly.
