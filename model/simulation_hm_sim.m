% Herding model - as modelled by Strombom et al
clear
clc
close all

%%
tic

% Parameters

no_shp = 14; % no.of sheep
box_length = 250; % initial box as considered in that model 
rad_rep_s = 2; % radius of repulsion from other sheep
K_atr = 10; % 10, 15, 25 (n = 14, 20, 30 resp)
k_atr = 8; % 8, 12, 20 k nearest neighbours to interact with
k_alg = 2; % 2, 2, 2
vs = 1; % speed of sheep
v_dog = 1.5; % speed of dog

h = 0.5; % relative strength of proceeding in the previous direction (inertia)
rho_a = 2; % relative strength of repulsion from other agents
rho_d = 1; % relative strength of repulsion from the dog
e = 0.3; % relative strength of angular noise
c = 1.5; % relative strength of attraction to the n nearest neighbours
alg_str = 1.3; % it was 1.3 for MS analysis

f_n = rad_rep_s*(no_shp)^(2/3);
pc = rad_rep_s;
pd = rad_rep_s;
rad_rep_dog = pd+1; % radius of repulsion from dog

no_it = 800;
n_iter = 500;

%%

parfor i = 1:no_it

    [pos_s_dat, pos_d_dat, vel_s_dat, vel_d_dat, spd_d_dat, collect_t, ...
        drive_t, force_slow_t] = herding_model_sim(no_shp, box_length, ...
        rad_rep_s, rad_rep_dog, K_atr, k_atr, k_alg, vs, v_dog, h, rho_a, rho_d, e, c, alg_str, ...
        f_n, pd, pc, n_iter)

    pos_s(:,:,:,i) = pos_s_dat;
    pos_d(:,:,i) = pos_d_dat;
    vel_s(:,:,:,i) = vel_s_dat;
    vel_d(:,:,i) = vel_d_dat;
    spd_d(:,i) = spd_d_dat;
    collect(:,i) = collect_t;
    drive(:,i) = drive_t;
    force_slow(:,i) = force_slow_t;

end

n_n = struct('pos_s', pos_s, 'pos_d', pos_d, 'vel_s', vel_s, 'vel_d', vel_d, 'no_it', no_it, ...
    'spd_d', spd_d, 'collect_t', collect, 'drive_t', drive, 'force_slow', force_slow, ...
    'no_shp', no_shp, 'box_length', box_length, 'rad_rep_s', rad_rep_s, 'rad_rep_dog', rad_rep_dog, ...
    'k_alg', k_alg, 'K_atr', K_atr, 'k_atr', k_atr, 'vs', vs, 'v_dog', v_dog, 'h', h, ...
    'rho_a', rho_a, 'rho_d', rho_d, 'e', e, 'alg_str', alg_str,...
    'c', c, 'f_n', f_n, 'pd', pd, 'pc', pc, 'n_iter', n_iter);
fname = strcat('hm_sim_', num2str(no_shp), '_rd_',num2str(round(rad_rep_dog)),'.mat');
save(fname, '-struct', 'n_n')

toc

disp('Simulation complete')