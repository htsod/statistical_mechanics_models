## Small-world Effect

The small-world effect is a mathematical model used to demonstrate the phenomenon of six degrees of separation. The key parameters of this model are the number of nodes $L$, the number of neighboring edges per node $k$, and the probability of creating shortcut connections $p$.

### Two regimes of behavior

By varying the value of $p$, we can significantly alter the network's topology. When $p$ is small, the number of shortcuts is near zero, and the network behaves like a regular graph. As $p$ increases, the graph becomes dominated by shortcuts, resulting in a fundamentally different structure.

![dense_and_sparse](/six_degree_separation/figures/merge.png)

We can use this model to explain our intuition behind the concept of six degrees of separation. To do so, we introduce a measure of the shortest distance between nodes.


### Average shortest Path length

The average shortest path length measures the average "closeness" of all node pairs in the graph. As we increase $p$ from 0.001 to 0.1, the number of shortcuts grows, and we track how this affects the average shortest path length.


![shortest_path](/six_degree_separation/figures/scaling_law.png)


The plot shows two distinct regimes of behavior as $p$ increases. This computational simulation supports our understanding of the six degrees of separation.

### Extension

For a more detailed analysis using the Renormalization method, refer to this [blog_post](https://htsod.github.io/posts/small_world/).




