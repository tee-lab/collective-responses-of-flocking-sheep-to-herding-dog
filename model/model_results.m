%% Run this code to produce Fig 6 of the main text

close all
clear
clc

%%

load('hm_n_14.mat')
fs = 30;
fname = 'Arial';

%% Calculating group polarisation (Fig 6c)

m = nan(n_iter, no_it);

for i = 1:no_it

    vel_temp = vel_s(:,:,:,i); % vel
    vel_temp = vel_temp./vecnorm(vel_temp, 2, 2); % norm vel
    m_temp = mean(vel_temp, 1); % mean vel
    m_temp = vecnorm(m_temp, 2, 2); % norm mean vel
    m(:,i) = squeeze(m_temp);
    
end

edges = 0:0.04:1;
m = m(:);
[prob_den_m, m_edges] = histcounts(m, edges, 'Normalization','pdf');
m_edges = m_edges(1:end-1) + (m_edges(2) - m_edges(1))/2;

model_fig_pol = figure('Position', [300 300 800 600]);

plot(m_edges, prob_den_m, '-' ,'Color', '#2365C4', 'LineWidth', 2.5)
set(gca, 'XLim', [0 1.02], 'XTick', 0:0.2:1, 'YLim', [0 7], 'YTick', 0:8, ...
    'FontName', fname, 'FontSize', fs, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Group polarisation', 'FontName', fname, 'FontSize', fs)
ylabel('PDF', 'FontName', fname, 'FontSize', fs)
exportgraphics(model_fig_pol, 'model_fig_pol.pdf', 'ContentType', 'vector')

%% Calculate group cohesion and distance of dog to barycenter

grp_coh = nan(n_iter, no_it);
dg_dist_bc = nan(n_iter, no_it);

for i = 1:no_it

    pos_sx = squeeze(pos_s(:,1,:,i));
    pos_sy = squeeze(pos_s(:,2,:,i));
    grp_cent_x = mean(pos_sx,1);
    grp_cent_y = mean(pos_sy,1);
    
    % calculate distance of sheep to barycenter along x and y
    pos_sx = pos_sx - grp_cent_x;
    pos_sy = pos_sy - grp_cent_y;

    pos_dx = pos_d(:,1,i);
    pos_dy = pos_d(:,2,i);
    
    % average distance of sheep to barycenter (cohesion) 
    dist_cent = sqrt(pos_sx.^2 + pos_sy.^2);
    dist_cent = mean(dist_cent,1);
    grp_coh(:,i) = dist_cent';

    % calculating distance of dog to barycenter
    pos_d_s_x = grp_cent_x - pos_dx';
    pos_d_s_y = grp_cent_y - pos_dy';
    dg_dist_bc(:,i) = sqrt(pos_d_s_x.^2 + pos_d_s_y.^2)';

end

grp_coh = grp_coh(:);
dg_dist_bc = dg_dist_bc(:);

gc_edges = min(grp_coh):0.2:max(grp_coh);
[gc_h, gc_edges] = histcounts(grp_coh, gc_edges, 'Normalization', 'pdf');
gc_edges = gc_edges(1:end-1) + (gc_edges(2) - gc_edges(1))/2;

dds_edges = min(dg_dist_bc):0.5:max(dg_dist_bc);
[dds_h, dds_edges] = histcounts(dg_dist_bc, dds_edges, 'Normalization', 'pdf');
dds_edges = dds_edges(1:end-1) + (dds_edges(2) - dds_edges(1))/2;

%% Calculating relative position of dog and sheep (Fig 6f)

r_sd_ci_all = nan(no_shp*(n_iter-1), no_it); %\psi_SD (angle between each sheep and dog)
r_dg_shp_ci_all = nan(no_shp*(n_iter-1), no_it); % \psi_DS (angle btw dog and each sheep)
r_bary_dg_ci_all = nan(n_iter-1,no_it); % \psi_BD
r_dg_bary_ci_all = nan(n_iter-1,no_it); % \psi_DB
r_phi_sd_all = nan(no_shp*(n_iter-1), no_it); % difference in orientation btw sheep and dog
r_phi_bary_dg_all = nan(n_iter-1,no_it); % \phi_BD

for i = 1:no_it

    pos_shp = pos_s(:,:,:,i); % sheep position
    vel_shp = vel_s(:,:,:,i); % sheep velovity
    phi_shp = squeeze(atan2(vel_shp(:,2,:), vel_shp(:,1,:))); % sheep orientation

    % dog position, velocity and orientation respectively.
    pos_dog = pos_d(:,:,i);
    pos_dog = pos_dog.';
    vel_dog = vel_d(:,:,i);
    phi_dog = atan2(vel_dog(:,2), vel_dog(:,1));

    tm = length(pos_shp);
    
    % storing all the above measures for a given simulation run
    r_sd_ci = zeros(no_shp, tm-1);
    r_dg_shp_ci = zeros(no_shp, tm-1);
    r_bary_dg_ci = zeros(1,tm-1);
    r_dg_bary_ci = zeros(1,tm-1);
    r_phi_sd = zeros(no_shp, tm-1);
    r_phi_bary_dg = zeros(1,tm-1);
    r_sd_sign = zeros(no_shp, tm-1);
    r_dg_shp_sign = zeros(no_shp, tm-1);
    vel_shp_sign = zeros(no_shp, tm-1);
    vel_dg_sign = zeros(no_shp, tm-1);

    for t = 2:tm

        r_sd_temp = repmat(pos_dog(:,t).', no_shp, 1) - pos_shp(:,:,t); % displacement vec between sheep and dog 
        r_bary_dg_temp = pos_dog(:,t).' - mean(pos_shp(:,:,t),1); % displacement vec between sheep bary and dog
        r_sd_temp = r_sd_temp./(sqrt(r_sd_temp(:,1).^2 + r_sd_temp(:,2).^2)); % normalising
        r_bary_dg_temp = r_bary_dg_temp./sqrt(r_bary_dg_temp(:,1).^2 + r_bary_dg_temp(:,2).^2); % normalising 
        r_dg_shp_temp = -r_sd_temp; % displacement vec between dog and sheep
        r_dg_bary_temp = -r_bary_dg_temp; % displacement vec between dog and sheep bary

        theta_shp = phi_shp(:,t); % orientation of sheep
        theta_dog = phi_dog(t); % orientation of dog
        vel_shp_temp = [cos(theta_shp) sin(theta_shp)]; % vel of sheep 
        vel_bary_temp = mean(vel_shp(:,:,t),1); % vel of bary
        vel_bary_temp = vel_bary_temp./sqrt(vel_bary_temp(:,1).^2 + vel_bary_temp(:,2).^2); % normalised vel of bary
        vel_dog_temp = [cos(theta_dog) sin(theta_dog)]; % vel of dog
        
        theta_bary_dg = dot(vel_bary_temp, vel_dog_temp, 2); 
        theta_bary_dg = acos(theta_bary_dg); % angle btw bary vel and vel of dog
        theta_bary_dg_sign = cross([vel_bary_temp 0], [vel_dog_temp 0],2); % calculating the sign of angle btw the 2
        theta_bary_dg = theta_bary_dg*sign(theta_bary_dg_sign(1,3));
        r_phi_bary_dg(1,t-1) = theta_bary_dg;  % angle btw sheep's barycentre and dog's orientation

        ci_bary_dg_temp = dot(r_bary_dg_temp, vel_bary_temp, 2); 
        ci_bary_dg_temp = acos(ci_bary_dg_temp); % angle btw bary and dog orientation
        ci_bary_dg_sign = cross([vel_bary_temp 0], [r_bary_dg_temp 0], 2);
        ci_bary_dg_temp = ci_bary_dg_temp*sign(ci_bary_dg_sign(1,3)); % sign of the angle
        r_bary_dg_ci(1,t-1) = ci_bary_dg_temp;  % where is dog from barycentre

        ci_dg_bary_temp = dot(r_dg_bary_temp, vel_dog_temp, 2); % angle btw dog and sheep from dog's orientation
        ci_dg_bary_temp = acos(ci_dg_bary_temp);
        ci_dg_bary_sign = cross([vel_dog_temp 0], [r_dg_bary_temp 0], 2);
        ci_dg_bary_temp = ci_dg_bary_temp*sign(ci_dg_bary_sign(1,3)); % sign og angle
        r_dg_bary_ci(1,t-1) = ci_dg_bary_temp;

        vel_dog_temp = repmat(vel_dog_temp, no_shp, 1);

        theta_sd = dot(vel_shp_temp, vel_dog_temp, 2);
        theta_sd = acos(theta_sd); % angle btw speed and dog

        ci_sd_temp = dot(r_sd_temp, vel_shp_temp, 2); % angle btw sheep and dog from sheeps orientation
        ci_dg_shp_temp = dot(r_dg_shp_temp, vel_dog_temp, 2); % angle btw dog and sheep from sheep's orientation
        ci_sd_temp = acos(ci_sd_temp);
        ci_dg_shp_temp = acos(ci_dg_shp_temp);
        r_sd_sign = [r_sd_temp zeros(no_shp,1)]; % calculating sign
        r_dg_shp_sign = [r_dg_shp_temp zeros(no_shp,1)];
        vel_shp_sign = [vel_shp_temp zeros(no_shp,1)];
        vel_dg_sign = [vel_dog_temp zeros(no_shp,1)];
        
        theta_sd_sign = cross(vel_shp_sign, vel_dg_sign, 2); % calculating sign of angle btw headings of sheep and dog
        theta_sd_sign = theta_sd_sign(:,3);
        theta_sd_sign = find(theta_sd_sign < 0);
        theta_sd(theta_sd_sign) = theta_sd(theta_sd_sign)*(-1);
        r_phi_sd(:,t-1) = theta_sd;

        ci_sd_sign = cross(vel_shp_sign, r_sd_sign, 2); % cal sign of angle between sheep's orientation and r_sd
        ci_dg_shp_sign = cross(vel_dg_sign, r_dg_shp_sign, 2);%  cal sign of angle between dog's orientation and r_ds
        ci_sd_sign = ci_sd_sign(:,3);
        ci_dg_shp_sign = ci_dg_shp_sign(:,3);
        ci_sd_sign = find(ci_sd_sign < 0);
        ci_dg_shp_sign = find(ci_dg_shp_sign < 0);

        ci_sd_temp(ci_sd_sign) = ci_sd_temp(ci_sd_sign)*(-1);
        ci_dg_shp_temp(ci_dg_shp_sign) = ci_dg_shp_temp(ci_dg_shp_sign)*(-1);

        r_sd_ci(:,t-1) = ci_sd_temp; % where is dog from sheep's prespective 
        r_dg_shp_ci(:,t-1) = ci_dg_shp_temp; % where is sheep from dog's perspective

    end

    r_sd_ci = r_sd_ci(:);

    r_dg_shp_ci = r_dg_shp_ci(:);

    r_phi_sd = r_phi_sd(:);

    r_sd_ci_all(:,i) = r_sd_ci(:);
    r_dg_shp_ci_all(:,i) = r_dg_shp_ci(:);
    r_bary_dg_ci_all(:,i) = r_bary_dg_ci(:);
    r_dg_bary_ci_all(:,i) = r_dg_bary_ci(:);
    r_phi_sd_all(:,i) = r_phi_sd(:);
    r_phi_bary_dg_all(:,i) = r_phi_bary_dg(:);

end

r_sd_ci = r_sd_ci_all(:);
r_dg_shp_ci = r_dg_shp_ci_all(:);
r_bary_dg_ci = r_bary_dg_ci_all(:);
r_dg_bary_ci = r_dg_bary_ci_all(:);
r_phi_sd = r_phi_sd_all(:);
r_phi_bary_dg = r_phi_bary_dg_all(:);

% converting angles to degrees from radians
r_sd_ci = (180*r_sd_ci)/pi;
r_dg_shp_ci = (180*r_dg_shp_ci)/pi;
r_bary_dg_ci = (180*r_bary_dg_ci)/pi;
r_dg_bary_ci = (180*r_dg_bary_ci)/pi;
r_phi_sd = (180*r_phi_sd)/pi;
r_phi_bary_dg = (180*r_phi_bary_dg)/pi;

edges = linspace(-180, 180, 50);

prob_den_bary_dg_ci = histcounts(r_bary_dg_ci, edges, 'Normalization', 'pdf');

prob_den_dg_bary_ci = histcounts(r_dg_bary_ci, edges, 'Normalization', 'pdf');

prob_den_phi_bary_dg = histcounts(r_phi_bary_dg, edges, 'Normalization', 'pdf');

edges = edges(1:(end-1)) + (edges(2) - edges(1))/2;

model_fig_psi_phi = figure('Position', [300 300 800 600]);
plot(edges, prob_den_bary_dg_ci, 'color', 'blue', 'LineWidth', 2)
hold on
plot(edges, prob_den_dg_bary_ci, 'color', 'red', 'LineWidth', 2)
hold on
plot(edges, prob_den_phi_bary_dg, '--', 'Color', '#0096FF', 'LineWidth', 2)

set(gca, 'XLim', [-180 180], 'XTick', -180:90:180, ...
    'FontName', fname, 'FontSize', fs, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Angle', 'FontName', fname, 'FontSize', fs)
ylabel('PDF', 'FontName', fname, 'FontSize', fs)

legend({'\psi_{\rm BD}', '\psi_{\rm DB}', '\phi_{\rm BD}'}, ...
    'Location', 'best', 'FontName', fname, 'FontSize', fs)
legend('boxoff')

exportgraphics(model_fig_psi_phi, 'model_fig_psi_phi.pdf', 'ContentType', 'vector')

%% Calculating group elongation

ry_rx = []; % group elongation
xd_final = []; % lateral displace of dog
yd_final = []; % rel pos of dog from barycenter along y axis along group velocity
ys_min_final = []; % y position of last sheep along group velocity

for i = 1:no_it

    pos_sheep = pos_s(:,:,:,i);
    vel_sheep = vel_s(:,:,:,i);
    pos_dog = pos_d(:,:,i)';

    grp_pos = mean(pos_sheep, 1);
    
    pos_shp_dg = pos_sheep;
    pos_shp_dg(15,:,:) = pos_dog;
    pos_shp_dg_cfr = pos_shp_dg - grp_pos;

    grp_vel = squeeze(mean(vel_sheep, 1));

    tm = size(vel_sheep,3); % length of the data set
    rx = nan(1,tm); % Length along the group velocity 
    ry = nan(1,tm); % length perpendicular to group velocity
    xd = nan(1,tm); % relative position of dog, xaxis
    yd = nan(1,tm); % relative position of dog, yaxis
    ys_min = nan(1,tm); % relative position of last sheep, yaxis.

    for t = 1:tm

        rot_ang = atan2(grp_vel(2,t), grp_vel(1,t)); % angle made by grp velocity along x-axis
        if rot_ang < 0
            rot_ang = rot_ang + 2*pi;
        end
        rot_mat = [cos(rot_ang) sin(rot_ang); -sin(rot_ang) cos(rot_ang)]; % rotation matrix
        pos_shp_cfr_temp = pos_shp_dg_cfr(:,:,t)'; % the vector should be of the for dimension X no.of particles
        new_pos = rot_mat*pos_shp_cfr_temp; % rotated coordinates
        rot_mat = [cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)];
        new_pos = rot_mat*new_pos;
        new_pos = new_pos';
        rx(t) = max(new_pos(1:no_shp,1)) - min(new_pos(1:no_shp,1)); % difference between max and min along x axis will the length along grp vel direction
        ry(t) = max(new_pos(1:no_shp,2)) - min(new_pos(1:no_shp,2)); % same for perpendicular to x grp vel direction
        xd(t) = new_pos(15,1);
        yd(t) = new_pos(15,2);
        ys_min(t) = min(new_pos(1:no_shp,2));    
    end

    ry_rx_temp = ry./rx;
    ry_rx = [ry_rx ry_rx_temp];
    xd_final = [xd_final xd];
    yd_final = [yd_final yd];
    ys_min_final = [ys_min_final ys_min];

end

edges = 0.1:0.2:7;

[ry_rx_hist, bin_edges] = histcounts(ry_rx, edges, 'Normalization', 'pdf');
    ryrx_bin_edges = bin_edges(1:end-1) + (bin_edges(2) - bin_edges(1))/2;

rel_yd = yd_final - ys_min_final;
rely_edges = -10:0.5:10;
[prob_rel_yd, rely_edges] = histcounts(rel_yd, rely_edges, 'Normalization','pdf');
rely_edges = rely_edges(1:end-1) + (rely_edges(2) - rely_edges(1))/2;

xd_rel_edges = -10:0.5:10;
[prob_xd_rel, xd_rel_edges] = histcounts(xd_final, xd_rel_edges, 'Normalization', 'pdf');
xd_rel_edges = xd_rel_edges(1:end-1) + (xd_rel_edges(2) - xd_rel_edges(1))/2;

%% Plotting group cohesion and group elongation (Fig 6 a)

model_fig_coh = figure('Position', [300 300 800 600]);
plot(gc_edges, gc_h, '-', 'Color', '#964B00', 'LineWidth', 2)

set(gca, 'XLim', [0 8], 'XTick', 0:1:8, 'YLim', [0 1.1], 'YTick', 0:0.2:1, ...
    'FontName', fname, 'FontSize', fs, 'LineWidth', 1)

xlabel('Cohesion (m)', 'FontName', fname, 'FontSize', fs);
ylabel('PDF', 'FontName', fname, 'FontSize', fs);

exportgraphics(model_fig_coh, 'model_fig_coh.pdf', 'ContentType', 'vector')

%% Plotting elongation (Fig 6 b)

model_fig_elong = figure('Position', [300 300 800 600]);
plot(ryrx_bin_edges, ry_rx_hist, '-', 'Color', '#A020F0', 'LineWidth', 2)

set(gca, 'XLim', [0 8], 'XTick', 0:1:8, 'YLim', [0 1.1], 'YTick', 0:0.2:1, ...
    'FontName', fname, 'FontSize', fs, 'LineWidth', 1)

xlabel('Elongation', 'FontName', fname, 'FontSize', fs);
ylabel('PDF', 'FontName', fname, 'FontSize', fs);

exportgraphics(model_fig_elong, 'model_fig_elong.pdf', 'ContentType', 'vector')

%% Plotting relative distance of dog to rear sheep (Fig 6 d)

model_fig_yrd_rbd = figure('Position', [300 300 800 600]); 
plot(rely_edges, prob_rel_yd, '-', 'Color', 'green', 'LineWidth', 2)

set(gca, 'XLim', [-10 10], 'XTick', -10:4:10, 'YLim', [0 .14], 'YTick', 0:0.02:.14, ...
    'FontName', fname, 'FontSize', fs, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Relative distance from dog to rare sheep', ...
    'FontName', fname, 'FontSize', fs)
ylabel('PDF', 'FontName', fname, 'FontSize', fs)
exportgraphics(model_fig_yrd_rbd, 'model_fig_yrd_rbd.pdf', 'ContentType', 'vector')

%% Plotting dog lateral movement (Fig 6e)

model_fig_xd = figure('Position', [300 300 800 600]); 
plot(xd_rel_edges, prob_xd_rel, '-', 'Color', '#FFD580' , 'LineWidth', 2)
set(gca, 'XLim', [-10 10], 'XTick', -10:4:10, 'YLim', [0 .085], 'YTick', 0:0.02:.1, ...
    'FontName', fname, 'FontSize', fs, 'LineWidth', 1)
xlabel('Dog lateral movement', 'FontName', fname, 'FontSize', fs)
ylabel('PDF', 'FontName', fname, 'FontSize', fs)
exportgraphics(model_fig_xd, 'model_fig_sd.pdf', 'ContentType', 'vector')
