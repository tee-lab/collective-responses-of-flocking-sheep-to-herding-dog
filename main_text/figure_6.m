%% Run this code to reproduce Fig.6a-b of the main text. 

close all
clear
clc

%% Loading data 

load('ext_drives_data.mat') 
no_evnts = length(events);

no_shp = no_ind - 2; % no.of sheep, i.e., no.of individuals - (dog + shepherd)
font_name = 'Arial';
font_size = 22;

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

mean_grp_coh = nan(1,size(gc_grp_spd_hs,2)); % mean cohesion for a given speed
std_err_grp_coh = nan(1,size(gc_grp_spd_hs,2)); % standard error

for i = 1:size(gc_grp_spd_hs,2)

    id_temp = bin_gs == i; % identify all the indices of speeds in a given bin
    gc_temp = grp_coh(id_temp); % find cohesion for corresponding indices
    mean_grp_coh(1,i) = mean(gc_temp);
    std_err_grp_coh(1,i) = (std(gc_temp)/sqrt(numel(gc_temp)));
    
end

%% calculate dog speed vs barycenter speed

[gs_ds_hs, spd_edges, ds_edges, bin_gs, bin_ds] = histcounts2(grp_spd, dog_spd', spd_edges, spd_edges);
 
mean_gs_ds = nan(1,size(gs_ds_hs,2)); % mean group speed for a given dog speed
std_err_gs_ds = nan(1,size(gs_ds_hs,2)); % std error

for i = 1:size(gs_ds_hs,2)

    id_temp = bin_ds == i; % identify all indicies for dog speed for a given bin
    gs_temp = grp_spd(id_temp); % corresponding group speed
    mean_gs_ds(1,i) = mean(gs_temp);
    std_err_gs_ds(1,i) = (std(gs_temp)/sqrt(numel(gs_temp)));
    
end

spd_edges = spd_edges(1:end-1) + (spd_edges(2) - spd_edges(1))/2;

%% Plotting Fig.6 of main text

fig_6 = figure('Position', [300 300 700 1200]);

subplot(2,1,1)

errorbar(spd_edges, mean_gs_ds, std_err_gs_ds, '-o', 'color', 'k',...
    'linewidth', 1.8, 'MarkerSize',10,'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k')
set(gca, 'XLim', [0 4], 'XTick', 0:.5:4, 'YLim', [0 2.5], 'YTick', 0:0.5:2.5, ...
    'FontSize', font_size, 'FontName', font_name, 'LineWidth', 1, 'Xcolor', 'k', 'YColor', 'k')
xl = xlabel('Dog speed (ms^{-1})', 'FontName', font_name, 'FontSize', font_size);
% xl.Position(2) = xl.Position(2) - abs(xl.Position(2));
yl = ylabel('Barycenter speed (ms^{-1})', 'FontName', font_name, 'FontSize', font_size);
% tl = title('a', 'FontSize', font_size+13, 'FontName', font_name, 'Color', 'k');
% tl.Position(1) = xl.Position(1) - 2.4;
% tl.Position(2) = tl.Position(2) - 0.1;

subplot(2,1,2)

errorbar(spd_edges, mean_grp_coh, std_err_grp_coh, '-o', 'color', "b",...
    'linewidth', 1.8, 'MarkerSize',10,'MarkerEdgeColor', "b", ...
    'MarkerFaceColor', "b")
set(gca, 'XLim', [0 4], 'XTick', 0:.5:4, 'YLim', [1 1.8], 'YTick', 1:0.2:1.8, ...
    'FontSize', font_size, 'FontName', font_name, 'LineWidth', 1, ...
    'Xcolor', 'k', 'YColor', 'k')
xl = xlabel('Barycenter speed (ms^{-1})', 'FontName', font_name, 'FontSize', font_size);
% xl.Position(2) = xl.Position(2) - abs(xl.Position(2)*0.07);
yl = ylabel('Cohesion (m)', 'FontName', font_name, 'FontSize', font_size);
% tl = title('b', 'FontSize', font_size+13, 'FontName', font_name, 'Color', 'k');
% tl.Position(1) = xl.Position(1) - 2.4;

exportgraphics(fig_6, 'fig_6.pdf', 'ContentType', 'vector')