![Ovis aries and border collie](main_text/P3160214.JPG)

# Collective responses of flocking sheep (*Ovis aries*) to a herding dog (border collie)

This repository contains the code for the manuscript:
Jadhav, V., Pasqua, R., Zanon, C., Roy, M., Tredan, G., Bon, R., Guttal, V. and Theraulaz, G., 2024. Collective responses of flocking sheep to a herding dog. bioRxiv, pp.2024-05. [bioRXiv preprint](https://www.biorxiv.org/content/10.1101/2024.05.24.595762v1.abstract)

The codes are tested to run on MATLAB Version: 23.2.0.2485118 (R2023b).

## Raw data

UWB tag data for all the trials are available in `/main_text/sheepR.dat`. The first, second, third, fourth and fifth columns in `sheepR.dat` are trial number, sheep ID, time (s), and position along the x (m) and y (m) axis, respectively. 

To convert the data in `sheepR.dat` into `.mat` format, run `/main_text/pre_ana.m` followed by `/main_text/pre_ini_ana.m`. This will give `sheep_all_dat.mat`. As our analysis focuses on the active phase, during which the dog was driving the sheep and staying behind the flock relative to its direction of movement, we select active phase trajectories from `sheep_all_dat.mat`. To do this, run `/main_text/drivs_dat.m`, which generates `drives_data.mat` 

Both data files `sheep_all_dat.mat` and `drives_data.mat` should be in both`\main_text` and `sm` folders to reproduce figures in those folders. We have already added these files to both the folder. 

## Code for Figure 1b

Run `\main_text\figure_1b` to reproduce representative trajectories of sheep, barycenter and dog shown in Figure 1b. 

## Codes for Figure 2

**Time series and probability density functions of the observables characterizing collective behavior of sheep and their reaction to the dog:** run `main_text\figure_2ae.m` to generate Figure 2a-e and run `\main_text\figure_2fj.m` to generate Figure 2f-j of the main text. 

## Code for Figure 3
**Relative angular positions and headings of the barycenter of the flock and the dog:** Run `main_text\figure_3.m` to generate Figure 3 of the main text. 

## Code for Figure 4
**Turning rates of the barycenter of the flock and the dog:** run `\main_text\figure_4.m` to generate Figure 4 of the main text. Ensure that `Violin.m` and `violinplot.m` (Bechtild, Bastian, 2016) are in path. 

## Codes for Figure 5
 **Hierarchical leader-follower relationships**
1. Run `\main_text\figure_5c` to generate leader-follower networks for all drives. Set variable `no_shp_dg = no_shp` to observe the leader-follower network without the dog and comment the above variable to observe the leader-follower network with the dog. This code generates Figure 5c and Supplementary Figure 5. 
2. Run `\main_text\figure_5de` to generate Figure 5 d-e of the main text and Supplementary Figures 6 and 7. To generate Fig 5e, the relationship between mean $d_i$ and indegree for the model, simulation data is already uploaded in the folder `main_text`. However, to run the simulations and calculate Pearson's correlation between mean $d_i$ and indegree, do the following:
	1. Run `\sm\simulation_hm.m`. Simulation data is stored as `hm_n_14.mat` for group size 14. 
	2. Load `hm_n_14.mat` in `\sm\model_network_ana.m` and run the code. 

## Codes for Figure 6
**Simulation results of the shepherding model**
To generate Figure 6 of the main text do the following:
1.  Run `\model\simulation_hm.m`. Simulation data is stored as `hm_n_14.mat` for group size 14. 
2. Load `hm_n_14.mat` in `\sm\model_results.m` and run the code.

## Modified herding model (Supplementary Section 3)

To simulate the modified model described in Supplementary Section 3, run `sm\simulation_hm_sim.m`. Simulation data is stored as `hm_sim_14_rd_3.mat` for $N = 14$ and radius of repulsion from the dog , $R_{\rm D} = 3$. For $N = 20$ and $30$, change model parameters as suggested in the code (see, Supplementary Fig.9 caption and Supplementary Table 1). To generate Supplementary Figure 9, run `\sm\model_network_ana.m` with data from appropriate $N$.   

## Supplementary Figures

Run codes in the folder `\sm` to generate Supplementary Figures. Ensure that `sheep_all_dat.mat` and `drives_data.mat` are in the path. 
1. `\sm\sm_fig_spd_coh.m` generates Supplementary Figure 2 
2. `\sm\sm_fig_cc.m` generates Supplementary Figures 3, 4, and 16.
3. `\sm\sm_fig_dij.m` generates Supplementary Figure 8 and 11. 
4. `\sm\all_trajectories.m` generates Supplementary Figure 14.

## References
1. Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.4559847


