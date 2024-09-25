## KAM Theorem

In an __ergodic system__, the system explores the entire energy surface over time. This concept is foundational in statistical mechanics because it allows for averaging over the energy surface to derive macroscopic observables in a system.

> Statistical mechanics focuses on large systems because it is well known that many system with a few interacting particles are definitely not ergodic.

For instance, the gravitational interaction between the Earth, Jupiter, and the Sun forms a non-ergodic three-body system. This system does not explore the full $18$-dimensional phase space (if it did, Earth's climate would be chaotic and unsuitable for life). To demonstrate this, a simulation of Earth's trajectory over 100 years was performed, plotting its orbit in the $x$-$y$ plane.

![100_years_trajectory](/kam_theorem/figures/one_hundred_years_orbit.png)



The Earth follows a stable orbit, both in Cartesian coordinates and in a rotating reference frame, confirming that the system remains non-ergodic under these conditions. However, if Jupiter's mass is artificially increased by a factor of 1000, Earth's orbit becomes highly irregular.

![heavy_jupiter](/kam_theorem/figures/heavy_jupiter.png)

### Quantifying Orbital Stability

This leads to the natural question: how does the irregularity of the orbit relate to the mass ratio of the three bodies? The first challenge is to quantify the stability of the orbit. Once we achieve that, we can explore the relationship between orbital stability and mass ratios.


To simplify the problem, we reduce the dimensionality from an $18$-dimensional phase space to a $2$-dimensional phase space, focusing on Earth's $R_x$ position and $V_y$ velocity. These reduced coordinates capture the system’s essential dynamics. Additionally, we record the reduced coordinates only when Earth crosses the line between the Sun and Jupiter, a method known as a Poincaré section. This reduction provides a clearer visualization of the orbital dynamics, often represented as tori in the reduced phase space. The first plot below shows the Poincaré section for the true mass of Jupiter.

![p_section(1)](/kam_theorem/figures/p_section(1).png)


As shown, the system forms a well-defined torus, indicative of regular, non-chaotic behavior. Next, we increase Jupiter's mass to 22,000 times that of Earth and observe a stark difference.

![p_section(2)](/kam_theorem/figures/p_section(2).png)


In this case, the system forms distinct, non-torus-like cycles. These deviations from a toroidal shape indicate irregular or chaotic motion. To further investigate, we plot Earth's orbit in a rotating frame of reference under these conditions.

![rot_22000](/kam_theorem/figures/rot_2200.png)

This phenomenon is called __mode-locking__. The system appears to be trapped in three distinct cycles, suggesting the presence of strange attractors that constrain the dynamics to specific modes based on the initial conditions.


### Sensitivity to Initial Conditions

We now explore how the system responds to different initial conditions. When varying these conditions, we find that the system behaves chaotically, with its final configuration highly sensitive to small changes in initial parameters—a hallmark of chaotic dynamics.


![p_section(3)](/kam_theorem/figures/p_section(3).png)

Each color in the plot represents a different initial condition. Notice the brown trajectory, which resembles a distorted torus, suggesting a circular orbit in Earth's trajectory. In topology, such distortions still preserve the essential properties of the original shape, even after twisting and stretching.

### The KAM Theorem and Topological Integrity

While we've made progress in visualizing the system's dynamics, the problem isn't fully solved. The KAM Theorem, developed by Kolmogorov, Arnold, and Moser, rigorously proves the stability of orbits in Hamiltonian systems under small perturbations. Their result shows that for sufficiently irrational winding numbers (such as the ratio of Jupiter’s orbital period to Earth’s), stability is preserved. More specifically, for sufficiently small planetary masses, the system's phase space contains a distorted torus, close to the unperturbed one, around which the planets spiral with the same winding number.

From the perspective of statistical mechanics, this three-body gravitational system presents a challenge to Boltzmann's ergodic hypothesis, which asserts that the system's states should equally sample the entire energy surface. Clearly, this system violates that assumption. This type of non-ergodicity can also be viewed through the lens of criticality. In small systems, the dynamical topology breaks continuous symmetries, creating asymmetrical potential wells. As more particles are added, the system becomes more ergodic, but in the limit of few particles, these asymmetries act as attractors.


### Extensions and Future Directions

One key insight from this problem is the utility of dimensional reduction and topological analysis in studying complex systems. By reducing the dimensionality and visualizing the topology, we can often simplify the problem and gain deeper understanding. Topological tools are becoming increasingly important in studying many-particle systems and could be further developed to justify coarse-graining and scaling procedures, which reduce complexity by decreasing dimensionality.