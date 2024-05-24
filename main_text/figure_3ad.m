%% Run this code to reproduce Fig.3a-d of the main text. 

close all
clear
clc

%% Load data

load('drives_data.mat')
load('sheep_all_dat.mat') % loading all sheep data --- pos, vel, etc
dr = 3; % We use drive 3 from event 9 to make this figure. 
no_shp = no_ind - 2; % no.of sheep, i.e., no.of individuals - (dog + shepherd)
font_name = 'Arial';
font_size = 30;
fig_pos = [300 300 800 400]; % position and size of figure

%% Calculate time series of group cohesion and distance of dog to barycenter

pos_x = squeeze(pos_ev_4(:,1,ev_st_4(dr):ev_et_4(dr))); % load x pos
pos_y = squeeze(pos_ev_4(:,2,ev_st_4(dr):ev_et_4(dr))); % load y pos

grp_pos_x = mean(pos_x(1:14,:),1); % group x pos
grp_pos_y = mean(pos_y(1:14,:),1); % group y pos

pos_s_x = pos_x(1:14,:) - grp_pos_x; % x_i - x_B
pos_s_y = pos_y(1:14,:) - grp_pos_y; % y_i - y_B

grp_coh = sqrt(pos_s_x.^2 + pos_s_y.^2); 
grp_coh = mean(grp_coh, 1);

dg_dist_bary_x = pos_x(15,:) - grp_pos_x; % x_D - x_B
dg_dist_bary_y = pos_y(15,:) - grp_pos_y; % y_D - y_B
dg_dist_bary = sqrt(dg_dist_bary_x.^2 + dg_dist_bary_y.^2); % distance of dog to barycenter

%% Calculate time series of group polarisation

% calculate mean and std dev of group polarisation of N (N = 14)
% non-interacting particles.

phi_temp = phi_ev_4(1:14,ev_st_4(dr):ev_et_4(dr)); % load heading angles (phi)
mx = mean(cos(phi_temp),1); % mx
my = mean(sin(phi_temp),1); % my
m = sqrt(mx.^2 + my.^2); % m

%% Calculate time series of elongation 

pos = pos_ev_4(:,:,ev_st_4(dr):ev_et_4(dr)); % load position
vel = vel_ev_4(:,:,ev_st_4(dr):ev_et_4(dr)); % load velocity

vel_sheep = vel(1:no_shp,:,:); % velocity of sheep
pos_sheep = pos(1:no_shp,:,:); % sheep position

grp_pos = mean(pos_sheep, 1); % grp center
grp_vel = squeeze(mean(vel_sheep, 1)); % grp velocity
pos_sheep_cfr = pos - grp_pos; % barycenter frame of reference

phi_B = atan2(grp_vel(2,:), grp_vel(1,:));

tm = size(vel_sheep,3); % length of the data set
rx = nan(1,tm); % Length along the group velocity
ry = nan(1,tm); % length perpendicular to group velocity
xd = nan(1,tm); % relative position of dog, xaxis
yd = nan(1,tm); % relative position of dog, yaxis
ys_min = nan(1,tm); % relative position of last sheep, yaxis.

for t = 1:tm

    rot_ang = atan2(grp_vel(2,t), grp_vel(1,t)); % angle made by grp velocity along x-axis
    rot_mat = [cos(rot_ang) sin(rot_ang); -sin(rot_ang) cos(rot_ang)]; % rotation matrix for group velocity along x axis
    pos_shp_cfr_temp = pos_sheep_cfr(:,:,t)'; % the vector should be of the form: dimension X no.of particles
    new_pos = rot_mat*pos_shp_cfr_temp; % rotated coordinates
    rot_mat = [cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)]; % rotate again by 90 degree so that group velocity is along y axis
    new_pos = rot_mat*new_pos;
    new_pos = new_pos';
    rx(t) = max(new_pos(1:no_shp,1)) - min(new_pos(1:no_shp,1)); % difference between max and min along x axis will the length perpendicular to grp vel direction
    ry(t) = max(new_pos(1:no_shp,2)) - min(new_pos(1:no_shp,2)); % same distance along grp vel direction
    xd(t) = new_pos(15,1); % distance of dog from barycenter along x-axis (\bar{x}_D in main text)
    yd(t) = new_pos(15,2); % distance of dog from barycenter along y-axis
    ys_min(t) = min(new_pos(1:no_shp,2)); % y distance of last sheep from barycenter. 

end

ry_rx = ry./rx; % elongated along group velocity if ry/rx > 1
y_rd = ys_min - yd; % relative distance from dog to rare sheep. +ve if dog is behind barycenter and -ve if dog is in front of the barycenter

%% Time series of speeds of sheep, barycenter and dog

shp_spd = squeeze(vecnorm(vel_sheep, 2, 2));
vel_bary = squeeze(mean(vel_sheep, 1));
spd_bary = vecnorm(vel_bary,2,1);
vel_dog = squeeze(vel(15,:,:));
spd_dog = vecnorm(vel_dog,2,1);
tm = (1:length(shp_spd))*dt;

%% Plotting Fig.3a of main text (ts of sheep, barycenter and dog speed respectively)

fig_3a = figure('Position', fig_pos);
t_stamps = ((51:50:length(pos_x))-1)*dt;

for i = 1:no_shp

    plot(tm, shp_spd(i,:), 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
    hold on

end

p1 = plot(tm, spd_bary, 'LineWidth', 2, 'Color', 'b');
hold on 
p2 = plot(tm, spd_dog, 'LineWidth', 2, 'Color', 'r');
hold on
xline(t_stamps, '--', 'LineWidth', 1)


set(gca, 'XLim', [0 max(tm)], 'XTick', 0:5:max(tm), 'YLim', [0 4], ...
    'YTick', 0:1:4, 'FontName', font_name, 'FontSize', font_size, ...
    'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Time (s)', 'FontName', font_name, 'FontSize', font_size)
ylabel('Speed','FontName', font_name, 'FontSize', font_size)

lg = [p1 p2];
legend(lg, 'v_B', 'v_D')
legend('boxoff')

%% Plotting Fig.3b of main text (ts of cohesion and elongation)

fig_3b = figure('Position', fig_pos);
p1 = plot(tm, grp_coh, '-', 'LineWidth', 2, 'Color', '#964B00');
hold on
p2 = plot(tm, ry_rx, '-.', 'LineWidth', 2, 'Color', '#A020F0');
hold on 
xline(t_stamps, '--', 'LineWidth', 1)

set(gca, 'XLim', [0 max(tm)], 'XTick', 0:5:max(tm), 'YLim', [0 2], ...
    'YTick', 0:0.5:2, 'FontName', font_name, 'FontSize', font_size, ...
    'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Time (s)', 'FontName', font_name, 'FontSize', font_size)
ylabel('Cohesion, E','FontName', font_name, 'FontSize', font_size)

lg = [p1 p2];
legend(lg, 'Cohesion', 'Elongation')
legend('boxoff')

%% Plotting Fig.3c of main text (ts of Polarisation)

fig_3c = figure('Position', fig_pos);
plot(tm, m, 'LineWidth', 2, 'Color', '#00008B')
hold on
xline(t_stamps, '--', 'LineWidth', 1)

set(gca, 'XLim', [0 max(tm)], 'XTick', 0:5:max(tm), 'YLim', [0.6 1], ...
    'YTick', 0.6:0.1:1, 'FontName', font_name, 'FontSize', font_size, ...
    'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Time (s)', 'FontName', font_name, 'FontSize', font_size)
ylabel('Polarization','FontName', font_name, 'FontSize', font_size)

%% Plotting Fig.3d of main text 
% (ts of relative distance of dog to rear sheep and distance of dog to barycenter)

fig_3d = figure('Position', fig_pos);
p1 = plot(tm, y_rd, '-.', 'LineWidth', 2, 'Color', 'green');
hold on
p2 = plot(tm, dg_dist_bary, '-', 'LineWidth', 2, 'Color', '#013220');
hold on 
xline(t_stamps, '--', 'LineWidth', 1)

set(gca, 'XLim', [0 max(tm)], 'XTick', 0:5:max(tm), 'YLim', [-0.5 8], ...
    'YTick', 0:2:8, 'FontName', font_name, 'FontSize', font_size, ...
    'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Time (s)', 'FontName', font_name, 'FontSize', font_size)
ylabel('$d_{\rm BD}$, $\bar{y}_{\rm RD}$', 'Interpreter', 'latex', ...
    'FontName', font_name, 'FontSize', font_size)

lg = [p1 p2];
legend(lg, 'Relative distance of dog to rare sheep', 'Distance of dog to barycenter', 'Location', ...
    'northwest')
legend('boxoff')

%% Plotting Fig.3e of main text (ts of Dog lateral movement and barycenter orientation)

fig_3e = figure('Position', fig_pos);
p1 = plot(tm, xd, '-', 'LineWidth', 2, 'Color', '#FFD580');
hold on
p2 = plot(tm, phi_B, '-.', 'LineWidth', 2, 'Color', '#FF7F7F');
hold on
xline(t_stamps, '--', 'LineWidth', 1)

set(gca, 'XLim', [0 max(tm)], 'XTick', 0:5:max(tm), 'YLim', [-4 5], ...
    'YTick', -4:2:5, 'FontName', font_name, 'FontSize', font_size, ...
    'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Time (s)', 'FontName', font_name, 'FontSize', font_size)
ylabel('$\bar{x}_{\rm D}$ (m), $\phi_{\rm B}$ (rad)', 'Interpreter', 'latex', ...
    'FontName', font_name, 'FontSize', font_size)

lg = [p1 p2];
legend(lg, 'Dog lateral movement', 'Baryceneter orientation', 'Location', ...
    'northwest')
legend('boxoff')