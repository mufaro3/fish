# Research Notes

Prof. Floryan left me the following two papers (links). My goal is to learn about the model he is using for the fish to be able to implement it as code, then simplify this algorithm based on the Barnes-Hut algorithm.

1. [Group cohesion and passive dynamics of a pair of inertial swimmers with three-dimensional hydrodynamic interactions.](https://arxiv.org/abs/2404.13152)

2. [Effects of symmetry and hydrodynamics on the cohesion of groups of swimmers.](https://doi.org/10.1088/1748-3190/ae0bd9)

## Notes

### 1. Dr. Floryan's Automaton-Dipole Fish Model

1. The fish are essentially modelled as two point-masses connected by a rod of length $\ell$.

2. There's a mathematical basis involving a "multipole expansion" (whatever that is) alongisde validations from other papers showing that using a source-sink pair separated by the length of the fish.

3. We orient them geometrically using symmetry to simplify them (essentially, we fix one of the fish and define the other fish relative to the first).

4. Their movement is (unsurprisingly) defined by a nonlinear first-order ODE.

5. The rest of the paper discusses the various resultant states that occur given different geometric configurations of the fish, as defined by three variables $\Delta \theta, \Delta x$, and $\Delta y$.

6. The fish move at a speed $$U = \frac{\sigma}{4\pi\ell^2}$$ where $\sigma \ge 0$ is the volumetric flow rate of the source and sink.

7. Each 'swimmer' $i$ has a center-of-mass position $\mathbf{x}_{c,i}$ and this has a separation between its source (head) at $\mathbf{x}_{f,i}$ and sink (tail) at $\mathbf{x}_{b,i}$ of $\ell$ along a direction denoted by the unit vector $\mathbf{n}_i$, related as:
$$\mathbf{x}_{f,i} = \mathbf{x}_{c,i} + \frac{1}{2} \ell \mathbf{n}_i$$
$$\mathbf{x}_{b,i} = \mathbf{x}_{c,i} - \frac{1}{2} \ell \mathbf{n}_i$$
this means that for a given time $t$, each swimmer's full orientation requires a minimum of five variables:
	1. The three-dimensional cartesian position of the fish, $\mathbf{x}_{c,i}$ (which is effectively three variables: $x_{c,i,x}, x_{c,i,y},$ and $x_{c,i,z}$). 
	2. The three-dimensional direction of the orientation unit vector $\mathbf{n}_i$, which can be reduced to two spherical dimensions $\rho_i$ and $\theta_i$ knowing that the vector must always remain a unit vector.

8. The dynamics of each fish is defined based on the aforementioned 5 free variables alongside the self-propelled speed $U$ to produce the self-propelled velocity of the source $\mathbf{v}_{f,i}$ as
$$\mathbf{v}_{f,i} = U \left( \mathbf{n}_i + \ell^2 \sum_{j\neq i}^N \left( \mathbf{r}_{f,i-f,j} - \mathbf{r}_{f,i-b,j} \right) \right)$$
and the sink $\mathbf{v}_{b,i}$ as
$$\mathbf{v}_{b,i} = U \left( \mathbf{n}_i + \ell^2 \sum_{j\neq i} (\mathbf{r}_{b,i-f,j} - \mathbf{r}_{b,i-b,j}) \right)$$
where
$$\mathbf{r}_{\alpha-\beta,i-j}=\frac{\mathbf{x}_{\alpha,i} - \mathbf{x}_{\beta,j}}{||\mathbf{x}_{\alpha,i} - \mathbf{x}_{\beta,j}||^3}$$

9. The translational ($\dot{\mathbf{x}}_{c,i}$) dynamics of the system is defined as
$$\dot{\mathbf{x}}_{c,i} = \frac{1}{2} (\mathbf{v}_{f,i} + \mathbf{v}_{b,i})$$

10. The rotational ($\dot{\mathbf{n}}_i$) dynamics of the system is defined as
$$\dot{\mathbf{n}}_i = \ell^{-1} ( \mathbf{v}_{f,i} - \mathbf{v}_{b,i} + 2\lambda_i \mathbf{n}_i )$$
where 
$$\lambda_i = -\frac{1}{2}[(\mathbf{v}_{f,i} - \mathbf{v}_{b,i}) \cdot \mathbf{n}_i]$$ 
is a Lagrange multiplier that ensures that $\ell$ is kept constant ($\dot{\ell} = 0$), and it originates 

11. With all of this in mind, we can concatenate all of the position $\mathbf{x}_{c,i}$ and orientation vectors $\mathbf{n}_i$ of the system to produce the state vector for the whole system, $\mathbf{X}$, and this vector has dimension $6N$ (or $5N$ using spherical definition for $\mathbf{n}$).

12. The system is evolved forward in time using fourth-order Runge-Kutta (i.e., computing $\dot{\mathbf{X}}$ for each time $t$ using a finite time-step $\delta t$ based on computing $\dot{\mathbf{x}}_{c,i}$ and $\dot{\mathbf{n}}_i$). 

### 2. The Barnes-Hut Algorithm

The fundamental problem of the project is that simulating large schools of fish is an $n$-body problem: the computational complexity (amount of time/compute required to calculate the motion) of the system scales on the order of $n^2$ ($\mathcal{O}(n^2)$), and systems of $n \ge 3$ are analytically unsolvable, requiring numerical simulation. This means that around $\log n \ge 2$ to $3$ ($n \approx 100$ to $1000$), these systems become too-complex to solve numerically as well (assuming we compute the differential evolution of each fish individually).

To solve this problem, we'll be employing the Barnes-Hut algorithm, which simplifies the system by clustering farther-out groups of fish within the school into a single conglomerate fish. This leads to increased (although marginal) error for a decrease in complexity to $\mathcal{O}(n \log n)$ 

Prof. Floryan left the following resources on the Barnes-Hut algorithm:

1. [The Barnes-Hut Approximation: Efficient computation of N-body forces](https://jheer.github.io/barnes-hut/)
2. [Fast Hierarchical Methods for the N-body Problem, Part 1](https://people.eecs.berkeley.edu/~demmel/cs267/lecture26/lecture26.html)
3. [The Barnes-Hut Algorithm](https://arborjs.org/docs/barnes-hut)
4. [A hierarchical O(N log N) force-calculation algorithm](https://www.nature.com/articles/324446a0)


Given a spatial distribution of fish as described in the state matrix $\mathbf{X}$ defined as
$$\mathbf{X} = (\mathbf{X}_1, \mathbf{X}_2, \dots)$$
where
$$\mathbf{X}_i = [\mathbf{x}_{c,i}, \mathbf{n}_i]^T,$$
the three-dimensional Barnes-Hut algorithm applied to fish would work as follows:

1. We determine an iteration order for the fish (in this case, we'll use the natural listing order described in $\mathbf{X}$), and choose the first fish, $\mathbf{X}_{1}$, to be the root of the Barnes-Hut Octree, $\mathcal{O}_{\mathbf{X}}$.
2. We iterate over $\mathbf{X}$ and insert each fish $\mathbf{X}$ according to spatial orientation, meaning that we subdivide the three-dimensional space recursively until all of the fish in the tree are simplified.
3. We recursively calculate the center of mass for each cell (both leaves and greater nodes).
4. For each fish in the tree, traverse the tree to compute the force on the fish through deciding whether or not to use the center of mass or individual points based on how far the fish yielding a contribution to the force are.
5. Use these forces to calculate dynamic change and evolve the state numerically, then repeat from step 1 for the next time-step.

There is also this paper for generating Octrees:

[Cornerstone: Octree Construction Algorithms for Scalable Particle Simulations](https://arxiv.org/abs/2307.06345)
