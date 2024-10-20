![Ovis aries and border collie](main_text/P3160214.JPG)

# Collective responses of flocking sheep (*Ovis aries*) to a herding dog (border collie)

This repository contains the code for the manuscript:
Collective responses of flocking sheep to a herding dog. Jadhav, V., Pasqua, R., Zano, C., Roy, M., Tredan G., Bon, R., Guttal, V., Theraulaz, G. (2024). The codes are tested to run on MATLAB Version: 23.2.0.2485118 (R2023b).

## Raw data

UWB tag data for all the trails is available in `/main_text/sheepR.dat`. The first, second, third, fourth and fifth columns in `sheepR.dat` are trail number, sheep ID, time, and position along the x and y axis, respectively. To convert the data in `sheepR.dat` into `.mat` format, run `/main_text/pre_ana.m` followed by `/main_text/pre_ini_ana.m`. This will give `sheep_all_dat.mat`. As our analysis focuses on the active phase, during which the dog was driving the sheep and staying behind the flock relative to its direction of movement, we select data from `sheep_all_dat.mat` only to select trajectories from the active phase. To do this, run `/main_text/drivs_dat.m`. `sheep_all_dat.mat` and `drives_data.mat` should be in `\main_text` and `sm` folders to reproduce figures in those folders. We have already added these files to those respective folders. 

## Code for Figure 1b

Run `\main_text\figure_1b` to reproduce representative trajectories of sheep, barycenter and dog shown in Figure 1b. 

## Codes for Figure 2

**Time series and probability density functions of the observables characterizing collective behavior of sheep and their reaction to the dog:** run `main_text\figure_2ae.m` to generate Figure 2a-e and run `\main_text\figure_2fj.m` to generate Figure 2f-j of the main text. 

## Code for Figure 3
**Relative angular positions and headings of the barycenter of the flock and the dog:** Run `main_text\figure_3.m` to generate Figure 3 of the main text. 

## Code for Figure 4
**Turning rates of the barycenter of the flock and the dog:** run `\main_text\figure_4.m` to generate Figure 4 of the main text. 

## Codes for Figure 5
 **Hierarchical leader-follower relationships**
1. Run `\main_text\figure_5c` to generate leader-follower networks for all drives. Set variable `no_shp_dg = no_shp` to observe the leader-follower network without the dog and comment the above variable to observe the leader-follower network with the dog. This code generates Figure 5c and Supplementary Figure 5. 
2. Run `\main_text\figure_5de` to generate Figure 5 d-e of the main text and Supplementary Figures 6 and 7. To generate Fig 5e, the relationship between mean $d_i$ and indegree for the model, model data is already uploaded in the folder `main_text`. However, to run the simulations and calculate Pearson's correlation between mean $d_i$ and indegree, do the following:
	1. Run `\sm\simulation_hm.m`. Simulation data is stored as `hm_n_14.mat` for group size 14. 
	2. Load `hm_n_14.mat` in `\sm\model_network_ana.m` and run the code. 

## Codes for Figure 6
**Simulation results of the shepherding model**
To generate Figure 6 of the main text do the following:
1.  Run `\sm\simulation_hm.m`. Simulation data is stored as `hm_n_14.mat` for group size 14. 
2. Load `hm_n_14.mat` in `\sm\model_results.m` and run the code. 


