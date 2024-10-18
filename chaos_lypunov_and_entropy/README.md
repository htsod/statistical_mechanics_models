## Chaos, Lyapunov, and Entropy Increase

### Solving the problem
This program mainly serves to demonstrate the character of chaos with Logistic Mapping. With $\mu = 0.9$, logistic mapping is known to demonstrate chaotic behavior. Meaning if we monitor two trajectories that were initially very close by, as we update the step, these trajectories eventually deviate exponentially. Here we numerically calculate the Lyapunov exponent to be $\mu \approx 0.2$. So, the deviation of the trajectories scale as:

$$ |\delta Z(t)| \approx e^{0.2 t} | \delta Z_{0} |  $$

The chaotic phenomena arises from this simple rules of logistic map is fasinating. It suggests that the final state of the chaotic system will decorrelate from its initial states, essentially destroy the information of the initial condition, resulting in a increase in entropy. We end this exercise with a quote from the book.


> Two configurations of classical atoms or billiard balls, with inital positions and velocities that are almost identical, will rapidly diverge as the collisions magnify small initial deviations in angle and velocity into large ones.

> It is this chaotic stretching, follding, and kneading of pahse space that is at the root of our explanation that entropy increases.

![separation](/chaos_lypunov_and_entropy/figures/separation.png)


### Possible extension to the problem

1. Detection for chaotic system

Suppose a system takes $n$ fixed parameter and have a evolution trajectory $x(t)$. We could in principle use the method in this problem to tell if the system behave chaotically base on the initial condition sensitivity properties of the chaotic system.


