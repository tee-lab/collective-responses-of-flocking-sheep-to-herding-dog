%% Run this code to reproduce Fig.3f-j of the main text. 

close all
clear
clc

%% Load data 

load('ext_drives_data.mat') % all drive data
no_evnts = length(events); % no.of events
no_shp = 14; % no.of sheep
font_name = 'Arial';
font_size = 30;
fig_pos = [300 300 800 400]; % position and size of figure

%% Calculating and plotting pdf of sheep, barycenter and dog speeds (Fig.3f)

dog_speed = [];
sheep_speed = [];
bary_speed = [];

for ev = 1:no_evnts

    vel = eval(strcat('vel_ev_',num2str(ev))); % load velocity
    vel_sheep = vel(1:14,:,2:end); % velocity of sheep
    vel_dog = squeeze(vel(15,:,2:end)); % velocity of dog
    
    bary_vel = mean(vel_sheep,1); % barycenter velocity
    bary_spd_temp = vecnorm(bary_vel,2,2); % barycenter speed
    bary_spd_temp = squeeze(bary_spd_temp);
    sheep_spd_temp = vecnorm(vel_sheep,2,2); % individual sheep speed
    sheep_spd_temp = sheep_spd_temp(:);
    dog_spd_temp = vecnorm(vel_dog,2,1); % dog speed
    
    % storing all data
    bary_speed = [bary_speed bary_spd_temp'];
    dog_speed = [dog_speed dog_spd_temp];
    sheep_speed = [sheep_speed sheep_spd_temp'];

end

spd_edges = 0:0.2:4; % edges to calculate pdfs

[shp_spd_hist, spd_edges] = histcounts(sheep_speed, spd_edges, 'Normalization', 'pdf');
[dg_spd_hist, spd_edges] = histcounts(dog_speed, spd_edges, 'Normalization', 'pdf');
[bary_spd_hist, spd_edges] = histcounts(bary_speed, spd_edges, 'Normalization', 'pdf');

spd_edges = spd_edges(1:end-1) + (spd_edges(2) - spd_edges(1))/2;

fig_3f = figure('Position', fig_pos);

plot(spd_edges, bary_spd_hist, '--', 'Color', 'magenta', 'LineWidth', 2)
hold on
plot(spd_edges, dg_spd_hist, '-', 'Color', '#de2d26', 'LineWidth', 2)
hold on
plot(spd_edges, shp_spd_hist, '-.', 'Color', [0.7 0.7 0.7], 'LineWidth', 2)

set(gca, 'XLim', [0 4], 'XTick', 0:.5:4, 'YLim', [0 1], ...
    'YTick', 0:0.2:1, 'FontName', font_name, 'FontSize', font_size, ...
    'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Speed (m/s)', 'FontName', font_name, 'FontSize', font_size)
ylabel('PDF','FontName', font_name, 'FontSize', font_size)

legend({'Barycenter speed', 'Dog speed', 'Sheep speed'}, 'FontName', font_name, ...
    'FontSize', font_size)
legend('boxoff')

mean_ds = mean(dog_speed)
std_ds = std(dog_speed)
se_ds = std_ds/sqrt(numel(dog_speed))

mean_ss = mean(sheep_speed)
std_ss = std(sheep_speed)
se_ss = std_ss/sqrt(numel(sheep_speed))

%% Calculating and plotting pdf of group polarisation (Fig.3h)

m_final = []; % group polarisation

for ev = 1:no_evnts

    phi_temp = eval(strcat('phi_ev_',num2str(ev))); % load heading angles (phi)

    mx = mean(cos(phi_temp(1:14,:)),1); % mx 
    mx = mx.';
    mx = mx(2:end); % not taking the 1st value because its always zero
    my = mean(sin(phi_temp(1:14,:)),1); % my
    my = my.';
    my = my(2:end); % not taking the 1st value because its always zero
    m = sqrt(mx.^2 + my.^2);
    m_final = [m_final; m];

end

% plotting pdf of group polarisation (m)
m_edges = 0:0.04:1;
[prob_den_m, m_edges] = histcounts(m_final, m_edges, 'Normalization','pdf');
m_edges = m_edges(1:end-1) + (m_edges(2) - m_edges(1))/2;

fig_3h = figure('Position', fig_pos);

plot(m_edges, prob_den_m, '-' ,'Color', '#00008B', 'LineWidth', 2)
set(gca, 'XLim', [0 1], 'XTick', 0:0.2:1, 'YLim', [0 7], 'YTick', 1:6, ...
    'FontSize', font_size, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Polarisation', 'FontName', font_name, 'FontSize', font_size);
ylabel('PDF', 'FontName', font_name, 'FontSize', font_size);

mode_m_id = find(prob_den_m == max(prob_den_m));
mode_m = m_edges(mode_m_id);
mean_m = mean(m_final)
std_m = std(m_final)
se_m = std_m/sqrt(numel(m_final))

%% Calculating pdf of group elongation  
 
ry_rx = []; % elongation
xd_final = []; % x position of dog wrt to barycenter
yd_final = []; % y position of dog wrt to barycenter
ys_min_final = []; % y position of last sheep

for ev = 1:no_evnts

    pos = eval(strcat('pos_ev_',num2str(ev))); % load position
    vel = eval(strcat('vel_ev_',num2str(ev))); % load velocity

    vel_sheep = vel(1:14,:,2:end); % velocity of sheep
    pos_sheep = pos(1:14,:,2:end); % sheep position

    grp_pos = mean(pos_sheep, 1); % barycenter position
    grp_vel = squeeze(mean(vel_sheep, 1)); % barycenter velocity
    pos_sheep_cfr = pos(:,:,2:end) - grp_pos; % all positions from barycenter reference frame
    
    tm = size(vel_sheep,3); % length of the data set
    rx = nan(1,tm); % Length along the group velocity 
    ry = nan(1,tm); % length perpendicular to group velocity
    xd = nan(1,tm); % relative position of dog, xaxis
    yd = nan(1,tm); % relative position of dog, yaxis
    ys_min = nan(1,tm); % relative position of last sheep, yaxis.

    for t = 1:tm

        rot_ang = atan2(grp_vel(2,t), grp_vel(1,t)); % angle made by grp velocity along x-axis
        rot_mat = [cos(rot_ang) sin(rot_ang); -sin(rot_ang) cos(rot_ang)]; % rotation matrix so that group velocity is along x axis
        pos_shp_cfr_temp = pos_sheep_cfr(:,:,t)'; % the vector should be of the form dimension X no.of particles
        new_pos = rot_mat*pos_shp_cfr_temp; % rotated coordinates
        rot_mat = [cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)]; % rotation matrix such that v_B is along yaxis
        new_pos = rot_mat*new_pos;
        new_pos = new_pos';
        rx(t) = max(new_pos(1:no_shp,1)) - min(new_pos(1:no_shp,1)); % difference between max and min along x axis (perpendiculat to v_B)
        ry(t) = max(new_pos(1:no_shp,2)) - min(new_pos(1:no_shp,2)); % difference btw max and min along yaxis (parallel to v_B)
        xd(t) = new_pos(15,1); % lateral position of dog
        yd(t) = new_pos(15,2); % \bar{y}_D
        ys_min(t) = min(new_pos(1:no_shp,2)); % position of the last sheep
    end

    ry_rx_temp = ry./rx; % elongation
    ry_rx = [ry_rx ry_rx_temp];
    xd_final = [xd_final xd]; % lateral position of dog
    yd_final = [yd_final yd]; % y position of dog wrt to barycenter
    ys_min_final = [ys_min_final ys_min]; % y position of last sheep
    
end

ryrx_bin_edges = 0.1:0.2:7; % elongation edges

[ry_rx_hist, ryrx_bin_edges] = histcounts(ry_rx, ryrx_bin_edges, 'Normalization', 'pdf');
ryrx_bin_edges = ryrx_bin_edges(1:end-1) + (ryrx_bin_edges(2) - ryrx_bin_edges(1))/2;

mode_ryrx_id = find(ry_rx_hist == max(ry_rx_hist));
mode_ryrx = ryrx_bin_edges(mode_ryrx_id)
mean_ryrx = mean(ry_rx)
median_ryrx = median(ry_rx)
std_ryrx = std(ry_rx)
se_ryrx = std_ryrx/sqrt(numel(ry_rx))

%% calculating pdf of group cohesion

grp_coh_final = []; % cohesion
dog_dist_bery_final = []; % distance from dog to barycenter

for ev = 1:no_evnts

    pos = eval(strcat('pos_ev_', num2str(ev))); % load position data
    pos_x = squeeze(pos(1:14,1,:)); % select positions of only sheep (x)
    pos_y = squeeze(pos(1:14,2,:)); % select positions of only sheep (y)
    grp_cent_x = mean(pos_x, 1); % calculating group centre (x)
    grp_cent_y = mean(pos_y, 1); % calculating group centre (y)

    pos_x = pos_x - grp_cent_x; % subtracting pos_x from grp_cent(x) (now everything in bary cent frame)
    pos_y = pos_y - grp_cent_y; % subtracting pos_x from grp_cent(y)

    grp_coh = sqrt(pos_x.^2 + pos_y.^2); % measuring the distance of each sheep from the group centre
    grp_coh = mean(grp_coh, 1); % avg dist from grp centre
    grp_coh_final = [grp_coh_final; grp_coh'];

    dog_dist_x = squeeze(pos(15,1,:)).' - grp_cent_x; % x distance of dog to barycenter
    dog_dist_y = squeeze(pos(15,2,:)).' - grp_cent_y; % y distance of dog to barycenter
    dog_dist_bery = sqrt(dog_dist_x.^2 + dog_dist_y.^2); % distance of dog to bary centre
    dog_dist_bery_final = [dog_dist_bery_final; dog_dist_bery'];

end

[prob_den_c, c_edges] = histcounts(grp_coh_final, 31, 'Normalization','pdf');
c_edges = c_edges(1:end-1) + (c_edges(2) - c_edges(1))/2;

[prob_den_db, db_edges] = histcounts(dog_dist_bery_final, 31, 'Normalization','pdf');
db_edges = db_edges(1:end-1) + (db_edges(2) - db_edges(1))/2;

mean_dB = mean(grp_coh_final)
std_dB = std(grp_coh_final)
se_dB = std_dB/sqrt(numel(grp_coh_final))

%% Calculating pdfs of relative distance between dog and barycenter

rel_yd = ys_min_final - yd_final; % relative distance of dog to rear sheep
[prob_rel_yd, rely_edges] = histcounts(rel_yd, 60, 'Normalization','pdf');
rely_edges = rely_edges(1:end-1) + (rely_edges(2) - rely_edges(1))/2;

[prob_xd_rel, xd_rel_edges] = histcounts(xd_final, 50, 'Normalization', 'pdf'); % dog lateral distance pdf
xd_rel_edges = xd_rel_edges(1:end-1) + (xd_rel_edges(2) - xd_rel_edges(1))/2;

%% Plotting pdf of cohesion and elongation Fig.3g

fig_3g = figure('Position', fig_pos);

plot(c_edges, prob_den_c, '-', 'LineWidth', 2, 'Color', '#964B00')
hold on
plot(ryrx_bin_edges, ry_rx_hist, '--', 'LineWidth', 2, 'Color', '#A020F0')

set(gca, 'XLim', [0 3.5], 'XTick', 0:0.5:3.5, 'YLim', [0 2.35], 'YTick', 0:0.5:2, ...
    'FontSize', font_size, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Cohesion (m), Elongation', 'FontName', font_name, 'FontSize', font_size);
ylabel('PDF', 'FontName', font_name, 'FontSize', font_size);

legend({'Cohesion', 'Elongation'}, 'FontName', font_name, 'FontSize', font_size)
legend('boxoff')

%% pdf of relative distance of dog to rare sheep and distance from dog to barycenter (Fig.3i)

fig_3i = figure('Position', fig_pos);

plot(rely_edges, prob_rel_yd, '-', 'LineWidth', 2, 'Color', 'green')
hold on
plot(db_edges, prob_den_db, '--', 'LineWidth', 2, 'Color', '#013220')

set(gca, 'XLim', [-4 12], 'XTick', -4:2:12, 'YLim', [0 0.4], 'YTick', 0:0.1:.4, ...
    'FontSize', font_size, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('$d_{\rm BD}$, $\bar{y}_{\rm RD}$ (m)', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('PDF', 'FontName', font_name, 'FontSize', font_size);

legend({'Relative distance from dog to barycenter', 'Distance from dog to baryceneter'}, ...
    'FontName', font_name, 'FontSize', font_size)
legend('boxoff')

%% pdf of dog lateral movement (Fig.3j)

fig_3j = figure('Position', fig_pos);

plot(xd_rel_edges, prob_xd_rel, 'LineWidth', 2, 'Color', '#FFD580')

set(gca, 'XLim', [-8 8], 'XTick', -8:2:8, 'YLim', [0 0.17], 'YTick', 0:0.05:.2, ...
    'FontSize', font_size, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('$\bar{x}_{\rm D}$ (m)', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('PDF', 'FontName', font_name, 'FontSize', font_size);

