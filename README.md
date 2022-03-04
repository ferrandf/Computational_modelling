# Computational_modelling

## First project: Active Carbon

### First assigment:
  
The convection by air of pollutant, with adsorption and desorption, in an Active Carbon (AC) filter can be modelled by the global problem (at canister or bulk scale) (for equations see report).

#### `main1LocalProblem.m`

1(a): The Matlab code `main1LocalProblem.m` solves the local problem for an AC grain assuming constant exterior concentration `c_ext = 1300 gr/m3`. We'll use Euler method for solving the system of ODEs. The Euler method is an explicit method which can approximate the values of the mean concentration (qb) and concentration on the boundary (qR) of the bullet using desired time step. We can expect order of convergence 1 (O(dt) where dt is the used timestep). Our equations and method doesn't have a space step to discuss its order of convergence, so there is no meaning to the question.

1(b) One first interesting study is trying to know how the material parameters (D_p or Intraparticular Diffusion and B (constant related to the Boundary conditions)) affect the solution (mean and boundary concentrations). To do so, we can plot the Loading (Figures 1+3k) and Unloading (Figures 2+3k) or both (Figures 3+3k) for different D_p and B.

We expect to see a faster adsorption/desorption of pollutant as the intraparticular diffusion increases. Both mean and boundary concentrations get affected.

Figure 1.1:

![dp08b135](https://user-images.githubusercontent.com/32911477/154801902-e60eaae1-c3bc-4f67-8eef-0e5a043ef2d8.png) 

Figure 1.2:

![dp07b135](https://user-images.githubusercontent.com/32911477/154801898-f1435220-101a-4c21-831a-1dff7ae40889.png) 

As we can see, the second figure shows a greater D_p and due to that, the adsorption is significantly faster both mean and boundary concentrations. The first figure clearly explains how the bullet is behaving: initially the concentration in the boundary increases rapidly but the mean one is slower. In 200s it seems both stabilize at one value (qb ≈ qR ≈ 0.3183). Approximately at this value the bullet would be considered to be full of pollutant.
At 200s, we impose the exterior concentration to be zero. At this point, the bullet starts the desorption of the pollutant, and as expected, the mean value of concentration inside the bullet is the one to decrease slower.

Let's now see what happens if we change the constant B. 

Figure 1.3:

![dp08b3](https://user-images.githubusercontent.com/32911477/154803429-05953407-52e8-497c-ae18-25259c2bb5b4.png)

Figure 1.4:

![dp08b05](https://user-images.githubusercontent.com/32911477/154803423-52506cc9-292c-4137-b35e-0a276d324613.png)


A curious thing happens: a greater B increases the concentration rapidly on the boundary of the bullet, but shows similar increase regarding the mean concentration. The behaviour for the desorption is similar (rapid decrease of boundary concentration and slower for the mean concentration).

A small B shows slower increase in both concentrations (for values of B <≈ 0.5 the bullet doesn't have time to fill itself of pollutant with 200s). 


#### `main2FD1Dcanister.m`

The Matlab code `main2FD1Dcanister.m` solves the couples problem, that is the global problem and the local problem, for a 1D case, with an upwind approximation for the convective term: 
```
∂c/∂x|i = (c_i - c_{i-1})/∆x + O(∆x) (first-order convergence)
``` 
The lines corresponding to the adsorption and desorption of AC grains are commented, thus the code initially solves just the convection-diffusion equation.

2(a) Firstly, we can check the stability condition C = |v|∆t/∆x ≤ 1 (Courant Number) for the convective term, with an upwind approximation. We can easily see that incrementing the Courant Number results in unstability (trying dt = 10 (Figure 6) we can see uncoherent results). We can also see that decreasing the Courant Number far from 1 (for example dt = 2 -> C = 0.50, figure 7) we get what seems to be a good result, but does not reflect the reality of the problem. Note that we are considering diffusion to be nearly zero (we can change the diffusion parameter to see different behaviour), so we would not expect to get diffusion. BUT WE SEE IT!! That "diffusion" is numerical error and does not represent reality.

Figure 1.5: Optimal C = 1 reached with dt = 4:

![dt4_canister](https://user-images.githubusercontent.com/32911477/154813172-1bba3220-a6ee-4864-a53f-8065bffcb2b5.png)

Figure 1.6: C = 2.5, dt = 10:

![dt10_canister](https://user-images.githubusercontent.com/32911477/154813192-44e1dee6-3768-496e-89cc-b456db2e71a9.png)

Figure 1.7: C = 0.50, dt = 2:

![dt2_canister](https://user-images.githubusercontent.com/32911477/154813201-322b7117-98c1-495e-afcc-a2bb06af94b4.png)

To sum up, choosing ∆t is extremely important, and the Courant Number for explicit methods is a great indicator to at least take a closer look to what the solution looks like and question ourselves if what we are seeing is due to numerical error or real behaviour.


2(b) Let's try to modify the original upwind approximation and consider a centered approximation for the convective term: 
```
∂c/∂x|i = (c_{i+1} - c_{i-1})/∆x^2 + O(∆x^2) (second-order convergence)
```
The file that does this is `main2_centered.m`. We can check that no time step would be small enough to correct the unstability that appears:

Figure 1.8:

![dt001](https://user-images.githubusercontent.com/32911477/155094465-6dadd900-50dd-4463-82c1-bccf9575bdcf.png)


2(c) Returning to the initial upwind approximation, we now can implement the effect of the AC grains (local problem). Check `main2_locglo.m` for the detailed code. Lets find out if we can use the same dt as before.
We see that we need a larger Tfinal = 30000 (why?? Because now we are considering the effect of the active carbon inside the canister, so the pollutant should move relatively slow) and a smaller dt. In fact dt = 0.1 is the first one that shows something. This could be due to the interaction between pollutant and active carbon.

We can actually estimate the velocity of the pollutant front looking at the plots:

Figure 9 (only global problem -> Tfinal = 100s):

![dp08dt4_glo](https://user-images.githubusercontent.com/32911477/154915316-2ca0dac3-fde1-4754-b348-ccbc9e41e0f3.png)


Figure 1.10 (local and global problem -> Tfinal = 30000s):

![dp08dt01_locglo](https://user-images.githubusercontent.com/32911477/154915188-2880264b-e8b9-49b6-a262-50ff2c9e9463.png)


We can see that the pollutant is approximately at x = 0.078 at Tfinal. This gives us a velocity of v = 2.6000e-06 cm/s

With that velocity we can estimate the time that takes to the pollutant to reach the end of the canister: t = 0.1/v ≈ 38462s

Plotting the problem for this amount of time should be sufficient to observe the front reaching x = 0.1: 


Figure 1.11:

![endcanister](https://user-images.githubusercontent.com/32911477/154917759-427510b9-4007-408a-bc25-c91cda02eca2.png)

2(d) Let's study the effect of the material and other constant parameters:

- q_m (initially at 0.4): Is the maximum concentration capacity inside the canister (at q = 0.4 the canister is full). Increasing the capacity of the canister should give us a slower transport of pollutant:

FIgure 1.12:

![qm08_locglo](https://user-images.githubusercontent.com/32911477/154919113-190d9902-f6f2-4449-96e0-b4721b3f091c.png)

We can clearly see that compared to figure 10, the front at Tfinal = 30000s is much more far away from the end of the canister when the q_m is increased (q_m = 0.8).

When we reduce the qm = 0.2, the end of the canister is already reached at T = 30000s:

Figure 1.13:

![qm02_locglo](https://user-images.githubusercontent.com/32911477/154919726-2c56612d-ba32-495b-87ba-a1418497925a.png)

- K (initially at 3.e-3): This constant is also used when we compute the Longmuir Isoterm (the grater the K the smaller becomes L(q)). We should understand the Longmuir Isoterm as the concentration inside the bullet. With this, we can expect that if we increase the K, then we have more space for pollutant, and so the canister would be full slower:

Figure 1.14 (K = 10.e-3):

![K10em3](https://user-images.githubusercontent.com/32911477/154921576-637a81ba-a5d5-4494-b59f-d3459f20472a.png)


- D_p (initially at 1.e-8): The D_p parameter changes the diffusion inside the bullet, but from a global point of view this doesn't affect the convection of pollutant through the canister. Changing the order of the parameter (D_p = 1.e-5) the plot seems to be the same:




                                                                           
- B (initially at 1.3514): In the local problem the B is 



2(e) It is clear that for the original values of the material parameters, the diffusion is not highly relevant in our solution. Its effect is included in the computation of `\sigmaconstant = ((1-epse)/epse)*3*Dp*B/R`. This constant is equal to 3.4514e-05 with the initial material parameters values. 
If we, for instance, change the value of D_p close to zero, (strage things happens, D_p = 1.e-10, dt = 0.1, Tfinal = 30000s)!!!!!!!

I believed that the diffusion doesn't effect the global problem (if I compute the main2 global problem without the local one, I don't see a significal change when I modify the value of D_p and B) but it has a small effect on the local one.


### Second assigment

#### `main_withoutNeumann.m`

In fluid dynamics a *potential flow* is a flow shuch that the velocity **v** can be expressed as: **v** = ∇\phi for some scalar function \phi, refered as velocity potential.


1(a) The potential flow is indeed irrotational: the curl of any gradient is equal to zero.

1(b) TODO

1(c) We will consider the domain: Ω = [0 10] x [-3 3], with impermeable material on the top and bottom sides, and some impermeable circular inclusions, and a difference of potential of 10 on the left and right sides: \phi(10) = \phi(0) + 10.


Figure 2.1: Boundary (axis [0 10 -3 3] with circle 1: x = 2.6, y = 0.5, r = 0.4; circle 2: x = 4, y = -1.4, r = 0.5; circle 3: x = 7.6, y = 1.5, r = 0.6)

![boundary_plot](https://user-images.githubusercontent.com/32911477/156813233-c681fb16-de35-444b-9ecc-4c913878998c.png)


Figure 2.2: Mesh (Using ez4u)

![mesh_dom2](https://user-images.githubusercontent.com/32911477/156812454-d794f83e-557a-4fd3-b601-fb5c14980062.png)

Figure 2.3: Mesh with dirichlet boundary points marked

![mesh_boundary](https://user-images.githubusercontent.com/32911477/156816517-a1b6f178-cfdb-4cca-8800-fae934980114.png)


1(c) TODO

1(d) Using the file `main_withoutNeumann.m`, previously changing the sourceTerm to zero and considering Dirichlet conditions on x = 0 and x = 10 (u_D = x), we get the following solution:

Figure 2.4: FEM solution

![fem_sol](https://user-images.githubusercontent.com/32911477/156816513-7371c586-6517-4575-af52-6bb7e32e0edb.png)

Figure 2.5: Analytical solution

![analytical_sol](https://user-images.githubusercontent.com/32911477/156816507-cf4e186d-a962-4f32-bb21-16972a6b1535.png)


1(e) If we represent the gradient using quiver in the nodal (2.6) and in the integration points (2.7) we obtain:

Figure 2.6:

![quiv_nodal](https://user-images.githubusercontent.com/32911477/156816523-e0b5ea19-50df-4a9e-8908-a7be18bc91eb.png)

Figure 2.7: 

![quiv_int](https://user-images.githubusercontent.com/32911477/156816520-ade7e30f-71ed-4a14-8147-d0d6f8bb8fab.png)

![quiv_int_zoom](https://user-images.githubusercontent.com/32911477/156818304-236f5b02-e7d3-4b4d-a205-ea6ab089ce30.png)


1(f) Plotting the streamlines:

Figure 2.8:

![streamlines](https://user-images.githubusercontent.com/32911477/156818900-a7436bf3-6226-45a1-a595-a2c24b3b6523.png)




