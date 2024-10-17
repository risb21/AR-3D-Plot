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
        <br>
        $$\Large p \Rightarrow \text{"Surface passes through cell"}$$
        <br>

        $$\Large \exists (u, v) \in cell,  \dfrac{f(u)}{||f(u)||} \neq \dfrac{f(v)}{||f(v)||} \Leftrightarrow p \mid u,v \in \mathbb{R}^3$$
        
    - For each edge with an intersection, "binary search" could be used to narrow the intersection point such that $f(x, y, z) \approx 0$

    - Get normals of each intersection point on the surface
        <br>

        $$\Large f'(x, y, z) = \begin{bmatrix}\dfrac{\partial f(x, y, z)}{\partial x} & \dfrac{\partial f(x, y, z)}{\partial y} & \dfrac{\partial f(x, y, z)}{\partial z} \end{bmatrix} = n_i$$
        
    - Differentiate using the limit definition 
        <br>

        $$\Large \dfrac{\partial f(x, y, z)}{\partial x} =  <!-- \underset{h \to 0}{\text{lim}} --> \dfrac{f(x + h, y, z) - f(x - h, y, z)}{2h}$$
        
    - Using the points of intersection and the normals at those points, a singular point within the cell is solved for, using the following formula
        <br>
        
        $$\Large \underset{x\in\mathbb{R^3}}{\text{min}}\space E[x] = \sum_{i = 1}^{k}(n_i \cdot (x - p_i))^2$$
        <br>

        $$\Large p_i \Rightarrow \text{i-th point of intersection on an edge}$$
        <br>

        $$\Large n_i \Rightarrow \text{normal at }p_i$$
        
    - This can be represented as matrix operations
        <br>

        $$\Large \underset{x \in \mathbb{R}}{\text{min}} \space  E[x] = (Ax - b)^T(Ax - b)$$
        <br> 

        $$\Large A =
        \begin{bmatrix}
            n_{1_{x}} & n_{1_{y}} & n_{1_{z}} \\
            n_{2_{x}} & n_{2_{y}} & n_{2_{z}} \\
            \vdots & \vdots & \vdots \\
            n_{k_{x}} & n_{k_{y}} & n_{k_{z}} 

        \end{bmatrix}_{k \times 3}$$
        <br>

        $$\Large b =
        \begin{bmatrix}
            n_{1} \cdot p_{1} \\
            n_{2} \cdot p_{2} \\
            \vdots \\
            n_{k} \cdot p_{k}
        \end{bmatrix}_{k \times 1}$$
        <br>
    
    - In the end, you must solve for $x$, which gives the location of the point with least error, corresponding to all normals on intersection points in the cell. 
        <br>

        $$\Large Ax = b$$
    
    - <u>**NOTE**</u> : If all normals on intersection points in the cell point in the same direction, there may be infinite solutions for $x$, even residing outside the cell.
    - To fix this, we can add a small bias towards the "Center of Mass" or average of intersection points within the cell.
        <br>

        $$\Large b = 2.5\times10^{-3}$$
        <br>

        $$\Large p_{CoM} = \underset{i = 1}{\overset{k}{\text{avg}}} \space p$$
        <br>

        $$\Large A =
        \begin{bmatrix}
            n_{1_{x}} & n_{1_{y}} & n_{1_{z}} \\
            n_{2_{x}} & n_{2_{y}} & n_{2_{z}} \\
            \vdots & \vdots & \vdots \\
            n_{k_{x}} & n_{k_{y}} & n_{k_{z}} \\
            bias & 0 & 0 \\
            0 & bias & 0 \\
            0 & 0 & bias
        \end{bmatrix}_{(k + 3) \times 3}$$
        <br>

        $$\Large p =
        \begin{bmatrix}
            p_{1} \\
            p_{2} \\
            \vdots \\
            p_{k} \\
            p_{CoM} \\
            p_{CoM} \\
            p_{CoM}
        \end{bmatrix}_{(k+3) \times 1}$$
        <br>



#### 2. Meshing
