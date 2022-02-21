# Computational_modelling

## First project: Active Carbon

The convection by air of pollutant, with adsorption and desorption, in an Active Carbon (AC) filter can be modelled by the global problem (at canister or bulk scale) (for equations see report).

### `main1LocalProblem.m`

The Matlab code `main1LocalProblem.m` solves the local problem for an AC grain assuming constant exterior concentration `c_ext = 1300 gr/m3`. We'll use Euler method for solving the system of ODEs. The Euler method is an explicit method which can approximate the values of the mean concentration (qb) and concentration on the boundary (qR) of the bullet using desired time step. We can expect order of convergence 1 (O(dt) where dt is the used timestep). Our equations and method doesn't have a space step to discuss its order of convergence, so there is no meaning to the question.

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


### `main2FD1Dcanister.m`

The Matlab code `main2FD1Dcanister.m` solves the couples problem, that is the global problem and the local problem, for a 1D case, with an upwind approximation for the convective term: 
```
∂c/∂x|i = (c_i - c_{i-1})/∆x + O(∆x) (first-order convergence)
``` 
The lines corresponding to the adsorption and desorption of AC grains are commented, thus the code initially solves just the convection-diffusion equation.

Firstly, we can check the stability condition C = |v|∆t/∆x ≤ 1 (Courant Number) for the convective term, with an upwind approximation. We can easily see that incrementing the Courant Number results in unstability (trying dt = 10 (Figure 6) we can see uncoherent results). We can also see that decreasing the Courant Number far from 1 (for example dt = 2 -> C = 0.50, figure 7) we get what seems to be a good result, but does not reflect the reality of the problem. Note that we are considering diffusion to be nearly zero (we can change the diffusion parameter to see different behaviour), so we would not expect to get diffusion. BUT WE SEE IT!! That "diffusion" is numerical error and does not represent reality.

Figure 5: Optimal C = 1 reached with dt = 4:

![dt4_canister](https://user-images.githubusercontent.com/32911477/154813172-1bba3220-a6ee-4864-a53f-8065bffcb2b5.png)

Figure 6: C = 2.5, dt = 10:

![dt10_canister](https://user-images.githubusercontent.com/32911477/154813192-44e1dee6-3768-496e-89cc-b456db2e71a9.png)

Figure 7: C = 0.50, dt = 2:

![dt2_canister](https://user-images.githubusercontent.com/32911477/154813201-322b7117-98c1-495e-afcc-a2bb06af94b4.png)

To sum up, choosing ∆t is extremely important, and the Courant Number for explicit methods is a great indicator to at least take a closer look to what the solution looks like and question ourselves if what we are seeing is due to numerical error or real behaviour.


Let's try to modify the original upwind approximation and consider a centered approximation for the convective term: 
```
∂c/∂x|i = (c_{i+1} - c_{i-1})/∆x^2 + O(∆x^2) (second-order convergence)
```
The file that does this is `main2_centered.m`. We can check that no time step would be small enough to correct the unstability that appears:

Figure 8:

![dt4_centered](https://user-images.githubusercontent.com/32911477/154813953-008948de-1fa8-4a29-b4e2-c381dc81506e.png)

Returning to the initial upwind approximation, we now can implement the effect of the AC grains (local problem). Check `main2_locglo.m` for the detailed code. Let's find out if we can use the same dt as before.
We see that we need a larger Tfinal = 30000 (why?? Because now we are considering the effect of the active carbon inside the canister, so the pollutant should move relatively slow) and a smaller dt. In fact dt = 0.1 is the first one that shows something. This could be due to the interaction between pollutant and active carbon.

We can actually estimate the velocity of the pollutant front looking at the plots:

Figure 9 (only global problem -> Tfinal = 100s):

![dp08dt4_glo](https://user-images.githubusercontent.com/32911477/154915316-2ca0dac3-fde1-4754-b348-ccbc9e41e0f3.png)


Figure 10 (local and global problem -> Tfinal = 30000s):

![dp08dt01_locglo](https://user-images.githubusercontent.com/32911477/154915188-2880264b-e8b9-49b6-a262-50ff2c9e9463.png)


We can see that the pollutant is approximately at x = 0.078 at Tfinal. This gives us a velocity of v = 2.6000e-06 cm/s

With that velocity we can estimate the time that takes to the pollutant to reach the end of the canister: t = 0.1/v ≈ 38462s

Plotting the problem for this amount of time should be sufficient to observe the front reaching x = 0.1: 


Figure 11:

![endcanister](https://user-images.githubusercontent.com/32911477/154917759-427510b9-4007-408a-bc25-c91cda02eca2.png)

Let's study the effect of the material and other constant parameters:

- q_m (initially at 0.4): Is the maximum concentration capacity inside the canister (at q = 0.4 the canister is full). Increasing the capacity of the canister should give us a slower transport of pollutant:

FIgure 12:

![qm08_locglo](https://user-images.githubusercontent.com/32911477/154919113-190d9902-f6f2-4449-96e0-b4721b3f091c.png)

We can clearly see that compared to figure 10, the front at Tfinal = 30000s is much more far away from the end of the canister when the q_m is increased (q_m = 0.8).

When we reduce the qm = 0.2, the end of the canister is already reached at T = 30000s:

Figure 13:

![qm02_locglo](https://user-images.githubusercontent.com/32911477/154919726-2c56612d-ba32-495b-87ba-a1418497925a.png)

- K (initially at 3.e-3): This constant is also used when we compute the Longmuir Isoterm (the grater the K the smaller becomes L(q)). We should understand the Longmuir Isoterm as the concentration inside the bullet. With this, we can expect that if we increase the K, then we have more space for pollutant, and so the canister would be full slower:

Figure 14 (K = 10.e-3):

![K10em3](https://user-images.githubusercontent.com/32911477/154921576-637a81ba-a5d5-4494-b59f-d3459f20472a.png)








