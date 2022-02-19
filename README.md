# Computational_modelling

## First project: Active Carbon

The convection by air of pollutant, with adsorption and desorption, in an Active Carbon (AC) filter can be modelled by the global problem (at canister or bulk scale) (for equations see report).

The Matlab code `main1LocalProblem.m` solves the local problem for an AC grain assuming constant exterior concentration `c_ext = 1300 gr/m3`. We'll use Euler method for solving the system of ODEs. The Euler method is an explicit method which can approximate the values of the mean concentration (qb) and concentration on the boundary (qR) of the bullet using desired time step. We can expect order of convergence 1 (O(dt) where dt is the used timestep).

One first interesting study is trying to know how the material parameters (D_p or Intraparticular Diffusion and B (constant related to the Boundary conditions)) affect the solution (mean and boundary concentrations). To do so, we can plot the Loading (Figures 1+3k) and Unloading (Figures 2+3k) or both (Figures 3+3k) for different D_p and B.

We expect to see a faster adsorption/desorption of pollutant as the intraparticular diffusion increases. Both mean and boundary concentrations get affected.

Figure 1:

![dp08b135](https://user-images.githubusercontent.com/32911477/154801902-e60eaae1-c3bc-4f67-8eef-0e5a043ef2d8.png) 

Figure 2:

![dp07b135](https://user-images.githubusercontent.com/32911477/154801898-f1435220-101a-4c21-831a-1dff7ae40889.png) 

As we can see, the second figure shows a greater D_p and due to that, the adsorption is significantly faster both mean and boundary concentrations. The first figure clearly explains how the bullet is behaving: initially the concentration in the boundary increases rapidly but the mean one is slower. In 200s it seems both stabilize at one value (qb ≈ qR ≈ 0.3183). Approximately at this value the bullet would be considered to be full of pollutant.
At 200s, we impose the exterior concentration to be zero. At this point, the bullet starts the desorption of the pollutant, and as expected, the mean value of concentration inside the bullet is the one to decrease slower.

Let's now see what happens if we change the constant B. 

Figure 3:

![dp08b3](https://user-images.githubusercontent.com/32911477/154803429-05953407-52e8-497c-ae18-25259c2bb5b4.png)

Figure 4:

![dp08b05](https://user-images.githubusercontent.com/32911477/154803423-52506cc9-292c-4137-b35e-0a276d324613.png)



A curious thing happens: a greater B increases the concentration rapidly on the boundary of the bullet, but shows similar increase regarding the mean concentration. The behaviour for the desorption is similar (rapid decrease of boundary concentration and slower for the mean concentration).

A small B shows slower increase in both concentrations (for values of B <≈ 0.5 the bullet doesn't have time to fill itself of pollutant with 200s). 
