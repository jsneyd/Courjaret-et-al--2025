# **2D Calcium Modelling (Read Me)**
** There are 3 groups of code that simulate different calcium tunnelling responses: **
1) Files ending with _Less_Tunnel: simulation that shows ineffective calcium tunnelling.
2) Files ending with _No_Tunnel: simulation that shows the absence of calcium tunnelling.
3) Files ending with _Tunnel: simulation that shows the presence of calcium tunnelling.

** There are 3 main files to generate and plot the time series of the calcium model: **
a) Geom_Mesh_ (generates a 2D mesh polygon of the model)
b) Sim_ (simulate the calcium model based on the polygon specified in Geom_Mesh_)
c) Plot_TS_ (generate time series plots which track the specified mesh points in the model based on the results produced from Sim_.

** Running the Simulation **
Step 1: Open one of the Geom_Mesh_ files depending on what to simulate (choose from Geom_Mesh_Less_Tunnel, Geom_Mesh_No_Tunnel, or Geom_Mesh_Tunnel).
Step 2: Run the Geom_Mesh code which saves an output file that contains the P, E, T matrices of the cytoplasm and ER domains, and also the interior boundaries info which the cytoplasm and ER domains share. The geometry and mesh size can be changed as desired, but do take note of the possible error that might arise from conflicting interior boundary points between the ER and cytoplasm. The error can be rectified by choosing a different mesh size for the domains.
Step 3: Open the corresponding Sim_ file.
Step 4: The simulation time and conditions (i.e. addition of CPA, external Ca2+, and Histamine) can be adjusted as desired. Change nothing if it is to simulate the figures in the paper.
Step 5: Run the Sim_ code which saves an output file that contains the Geom_Mesh info and the calcium concentration at each mesh points throughout the simulation.
Step 6: Open the corresponding Plot_TS_ file.
Step 7: The tracked mesh points can be changed to reflect the calcium concentration to be tracked at a specific point in the polygon. Please refer to the geometric mesh figure generated in Geom_Mesh_ to identify the coordinates of the mesh points.
Step 8: Run the Plot_TS_ code which produces 2 time series plots tracking the changes of the calcium concentration simulated in the Sim_ code in the cytoplasm and in the ER at the specified points, respectively.

** Movie code (plot_tunnel_movie.m) **
This code is similar to the Plot_TS_ code except that it generates 2 movies, each tracking the calcium concentration in the cytoplasm and in the ER domains across time. The code utilises the same output file produced from the Sim_ code.

