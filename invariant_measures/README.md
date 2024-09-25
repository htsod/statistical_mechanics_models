## Invariant Measure

Here, we will be studying the evolution of non-hamiltonian system, which is again the logistic mapping. One reason this is a good math model to study is because by twigging the system parameter $\mu$, the evolution of the system displays different kind of dynamical behavior.

### Solving the problem

- Attracting fixed point

For $\mu < \frac{1}{4}$, the trajector has a fixed point at $x^{\ast}=0$

![fixed_point_zero](/invariant_measures/attracting_fixed_point(1).png)

If $1/4 <\mu < 3/4$, there is a non-zero fixed point.

![non-zero_fixed_point](/invariant_measures/attracting_fixed_point(2).png)


For $3/4 < \mu < (1+\sqrt{6})/4$, the trajectory has a stable, period-two cycle.

![period_two_cycle](/invariant_measures/period_two_cycle.png)

If $\mu = 1$, the trajectory explore all the configuration with unequal probability. That would suggest that the dynamics is ergodic but does not follow from Liouville's theorem.

![ergodic_dynamics](/invariant_measures/ergodic_dynamics.png)

Since each configuration has different weight, we would like to know the weight density, or calling it the invariant density.

![invariant_measure(1)](/invariant_measures/invariant_density(1).png)

It does cover the whole domain with weights tilt towards the ends. Things become tricky when $\mu$ is slightly less than $1$, cusps in the invariant density would appear.

![invariant_density(2)](/invariant_measures/invariant_density_mu(2).png)

![compare_cusps](/invariant_measures/compare_cusps.png)

Next, we will be studying how the invariant density varies as $\mu$ slightly increases from 0.8 to 1. Initially, we would observe a double peaks in the density suggesting there is a period two cycle. As $\mu$ increases, the invariant density doubles in its periodicity and eventually becomes noise. This is also known as the onset to chaos, where the final configuration depends sensitively to the initial configuration.

![bifurcation](/invariant_measures/variation_invariant_measure.png)

### Possible extension

The concept of the invariant measure is very useful when we studying a system with known configuration but with a unpredictable trajectories like that of the logistic mapping. It gives a qualitative description to a complicated system. We can thus design a numerical program that take in a function that descripts a particular system and output the corresponding invariant measure.


