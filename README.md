# AR-3D-Plot
AR application that recognizes 3D formulae and visualizes them in an AR environment

## Dual Contouring 
An algorithm to approximate implicit surfaces. This is used in this project to generate the 3D plots of the surfaces described by 3D formulae.

### References
1. Paper on [Dual Contouring on Hermite Data](https://www.cs.rice.edu/~jwarren/papers/dualcontour.pdf)
2. Blog describing [point generation in Dual Contouring](https://www.boristhebrave.com/2018/04/15/dual-contouring-tutorial/)
3. Blog describing [meshing points together after point generation](https://bonsairobo.medium.com/smooth-voxel-mapping-a-technical-deep-dive-on-real-time-surface-nets-and-texturing-ef06d0f8ca14#a230)

### Explanation
Dual Contouring comprises of 2 main steps:
1. Point Generation 
2. Meshing of points

#### 1. Point Generation 
- Sample points int 3D space at fixed and regular intervals, according to some function of all 3 axes, $f(x, y, z)$. The point will be on the surface if $f(x, y, z) = 0$ .
- Now traverse the samples in the form of 3D cubes/cells. Each cell considers 8 points at a time.
    - Each edge in the cell is checked for intersections, where an intersection has occurred if the sample al each end of the edge has a different sign.

    $$p \Rightarrow \text{"Surface passes through cell"}$$
    $$\exists (u, v) \in cell \rightarrow \dfrac{f(u)}{||f(u)||} \neq \dfrac{f(v)}{||f(v)||} \mid u,v \in \mathbb{R}^3 \Rightarrow p$$
    - For each edge with an intersection, "binary search" could be used to narrow the intersection point such that $f(x, y, z) \approx 0$
    - Get normals of each intersection point on the surface
    $$f'(x, y, z) = \left[\dfrac{\partial f(x, y, z)}{\partial x} \dfrac{\partial f(x, y, z)}{\partial y} \dfrac{\partial f(x, y, z)}{\partial z} \right] = n_i$$
    - Differentiate using the limit definition 
    <br>
    $$\dfrac{\partial f(x, y, z)}{\partial x} = \lim_{h\to 0} \dfrac{f(x + h, y, z) - f(x - h, y, z)}{2h}$$
    - Using the points of intersection and the normals at those points, a singular point within the cell is solved for, using the following formula

    $$\underset{x\in\mathbb{R^3}}{\text{min}}\space E[x] = \sum_{i = 0}^{k}(n_i \cdot (x - p_i))^2\newline$$

    $$p_i \Rightarrow \text{i-th point of intersection on an edge} \newline$$

    $$n_i \Rightarrow \text{normal at }p_i \newline$$
    - This can be represented as matrix operations

#### 2. Meshing
