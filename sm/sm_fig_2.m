%% Run this code to reproduce Supplementary Fig.2a-b.

close all
clear
clc

%% Loading data 

load('drives_data.mat') 
no_evnts = length(events);

no_shp = no_ind - 2; % no.of sheep, i.e., no.of individuals - (dog + shepherd)
font_name = 'Arial';
font_size = 22;

fig_2_sm = figure('Position', [300 300 1400 450]);

%% combining all data
    
phi_all_evts = []; % heading direction
vel_all_evts = []; % velocity
pos_all_evts = []; % position

for ev = 1:no_evnts

    phi_temp = eval(strcat('phi_ev_',num2str(ev))); % load heading angles (phi)
    vel_temp = eval(strcat('vel_ev_',num2str(ev))); % load velocity
    pos_temp = eval(strcat('pos_ev_', num2str(ev))); % load position data
    phi_all_evts = [phi_all_evts phi_temp(:,2:end)];
    vel_all_evts = cat(3, vel_all_evts, vel_temp(:,:,2:end)); % concatenate to the 3rd dimension
    pos_all_evts = cat(3, pos_all_evts, pos_temp(:,:,2:end));

end

%% Calculating barycenter speed vs cohesion

vel_sheep = vel_all_evts(1:14,:,:); % velocity of sheep
vel_dog = squeeze(vel_all_evts(15,:,:)); % velocity of dog

dog_spd = sqrt(vel_dog(1,:).^2 + vel_dog(2,:).^2); % dog speed
dog_spd = dog_spd(2:end); % considering only valid spds.

vel_grp = mean(vel_sheep,1); % barycenter velocity 
grp_spd = vecnorm(vel_grp,2,2); % barycenter speed
grp_spd = squeeze(grp_spd);
grp_spd = grp_spd(2:end);

spd_edges = 0:0.2:4; % speed bins

pos_x = squeeze(pos_all_evts(1:14,1,:)); % select positions of only sheep (x)
pos_y = squeeze(pos_all_evts(1:14,2,:)); % select positions of only sheep (y)
grp_cent_x = mean(pos_x, 1); % calculating group centre (x)
grp_cent_y = mean(pos_y, 1); % calculating group centre (y)

pos_x = pos_x - grp_cent_x; % subtracting pos_x from grp_cent(x) (now everything in bary cent frame)
pos_y = pos_y - grp_cent_y; % subtracting pos_x from grp_cent(y)

grp_coh = sqrt(pos_x.^2 + pos_y.^2); % measuring the distance of each sheep from the group centre
grp_coh = mean(grp_coh, 1); % avg dist from grp centre
grp_coh = grp_coh(2:end)';

% calculate pdf in 2d
gc_edges = 0.7:0.4:3.5;
[gc_grp_spd_hs, gc_edges, ~, bin_gc, bin_gs] = histcounts2(grp_coh, grp_spd, gc_edges, spd_edges);

spd_edges = spd_edges(1:end-1) + (spd_edges(2) - spd_edges(1))/2;

mean_grp_coh = nan(1,size(gc_grp_spd_hs,2)); % mean cohesion for a given speed
std_err_grp_coh = nan(1,size(gc_grp_spd_hs,2)); % standard error
subplot(1,2,2)

for i = 1:size(gc_grp_spd_hs,2)

    id_temp = bin_gs == i; % identify all the indices of speeds in a given bin
    gc_temp = grp_coh(id_temp); % find cohesion for corresponding indices
    mean_grp_coh(1,i) = mean(gc_temp);
    std_err_grp_coh(1,i) = (std(gc_temp)/sqrt(numel(gc_temp)));
    % gc_bp(1:length(gc_temp),i) = gc_temp;

    plot(spd_edges(i)*ones(length(gc_temp),1), gc_temp, 'o', 'Color', "#B9DDF1", ...
        'MarkerSize', 2, 'MarkerFaceColor', '#B9DDF1')
    hold on
    
end

errorbar(spd_edges, mean_grp_coh, std_err_grp_coh, '-o', 'color', "#3182bd",...
    'linewidth', 2, 'MarkerSize',8,'MarkerEdgeColor', "#3182bd",'MarkerFaceColor', "#3182bd")
set(gca, 'XLim', [0 4], 'XTick', 0:.5:4, 'YLim', [min(grp_coh) max(grp_coh)], ...
    'YTick', 0.5:0.5:max(grp_coh), 'FontSize', font_size, 'FontName', font_name, ...
    'LineWidth', 1, 'Xcolor', 'k', 'YColor', 'k')
xlabel('Barycenter speed (ms^{-1})', 'FontName', font_name, 'FontSize', font_size);
ylabel('Cohesion (m)', 'FontName', font_name, 'FontSize', font_size);

%% calculate dog speed vs barycenter speed

spd_edges = 0:0.2:4; % speed bins

[gs_ds_hs, spd_edges, ds_edges, bin_gs, bin_ds] = histcounts2(grp_spd, dog_spd', spd_edges, spd_edges);
spd_edges = spd_edges(1:end-1) + (spd_edges(2) - spd_edges(1))/2;
 
mean_gs_ds = nan(1,size(gs_ds_hs,2)); % mean group speed for a given dog speed
std_err_gs_ds = nan(1,size(gs_ds_hs,2)); % std error

subplot(1,2,1)

for i = 1:size(gs_ds_hs,2)

    id_temp = bin_ds == i; % identify all indicies for dog speed for a given bin
    gs_temp = grp_spd(id_temp); % corresponding group speed
    mean_gs_ds(1,i) = mean(gs_temp);
    std_err_gs_ds(1,i) = (std(gs_temp)/sqrt(numel(gs_temp)));
    plot(spd_edges(i)*ones(length(gs_temp),1), gs_temp, 'o', 'Color', '#C3C3C3', ...
        'MarkerSize', 2, 'MarkerFaceColor', '#C3C3C3')
    hold on
    
end

errorbar(spd_edges, mean_gs_ds, std_err_gs_ds, '-o', 'color', '#2B2B2B',...
    'linewidth', 2, 'MarkerSize',10,'MarkerEdgeColor', '#2B2B2B', ...
    'MarkerFaceColor', '#2B2B2B')
set(gca, 'XLim', [0 4], 'XTick', 0:.5:4, 'YLim', [0 max(grp_spd)], 'YTick', 0:1:max(grp_spd), ...
    'FontSize', font_size, 'FontName', font_name, 'LineWidth', 1, 'Xcolor', 'k', 'YColor', 'k')
xlabel('Dog speed (ms^{-1})', 'FontName', font_name, 'FontSize', font_size);
ylabel('Barycenter speed (ms^{-1})', 'FontName', font_name, 'FontSize', font_size);

exportgraphics(fig_2_sm, 'sm_fig_2.pdf', 'ContentType', 'vector')
% exportgraphics(fig_2_sm, 'sm_fig_2.jpeg')