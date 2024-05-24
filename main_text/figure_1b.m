%% Run this code to reproduce Fig.1b of the main text. 

close all
clear
clc

%% Load data

load('drives_data.mat')
load('sheep_all_dat.mat') % loading all sheep data --- pos, vel, etc

%% We use drive 3 from event 9 to make this figure. 

dr = 3;
pos_x = squeeze(pos_ev_4(:,1,ev_st_4(dr):ev_et_4(dr))); % load x position
pos_y = squeeze(pos_ev_4(:,2,ev_st_4(dr):ev_et_4(dr))); % load y position

vel_shp = vel_ev_4(:,:,ev_st_4(dr):ev_et_4(dr)); % load velocity
vel_shp_x = squeeze(vel_ev_4(:,1,ev_st_4(dr):ev_et_4(dr)));
vel_shp_y = squeeze(vel_ev_4(:,2,ev_st_4(dr):ev_et_4(dr)));

vel_shp_norm = squeeze(vecnorm(vel_shp,2,2)); % normalise the velocity
vel_shp_x = vel_shp_x./vel_shp_norm;
vel_shp_y = vel_shp_y./vel_shp_norm;

fig_1b = figure('Position', [300 300 800 600]);

% plotting trajectories of sheep

for i = 1:14
    
    if i == 1
       h1 = plot(pos_x(i,:), pos_y(i,:), '-o', 'Color', '#56B4E9', 'LineWidth', ...
           1, 'MarkerIndices', 1:50:length(pos_x), 'MarkerSize', 10, 'MarkerFaceColor', ...
           '#56B4E9', 'MarkerEdgeColor', '#56B4E9');
        hold on
    else
        plot(pos_x(i,:), pos_y(i,:), '-o', 'Color', '#56B4E9', 'LineWidth', ...
           1, 'MarkerIndices', 1:50:length(pos_x), 'MarkerSize', 10, 'MarkerFaceColor', ...
           '#56B4E9', 'MarkerEdgeColor', '#56B4E9')
        hold on
    end

end

hold on

% plotting trajectory of barycenter

h2_col = '#004f7c';

h2 = quiver(mean(pos_x(1:14,1:50:end),1), mean(pos_y(1:14,1:50:end),1), ...
    mean(vel_shp_x(:,1:50:end),1), mean(vel_shp_y(:,1:50:end),1), 0.17, 'color', h2_col, ...
    'Marker', '.', 'MarkerSize', 38, 'LineWidth', 2, 'MarkerEdgeColor', h2_col, ...
    'MarkerFaceColor', h2_col);

hold on

% plotting trajectory of dog

h3 = plot(pos_x(15,:), pos_y(15,:), '-o', 'Color', 'red', 'LineWidth', ...
           1, 'MarkerIndices', 1:50:length(pos_x), 'MarkerSize', 10, 'MarkerFaceColor', ...
           'red', 'MarkerEdgeColor', 'red');

hold on 

% plotting trajectory of shepherd

% h4 = plot(pos_x(16,:), pos_y(16,:), '-o', 'Color', '#009E73', 'LineWidth', ...
%            1, 'MarkerIndices', 1:50:length(pos_x), 'MarkerSize', 10, 'MarkerFaceColor', ...
%            '#009E73', 'MarkerEdgeColor', '#009E73');

lg_plt = [h1 h2 h3];

font_name = 'Arial';
font_size = 30;

legend(lg_plt, 'Sheep', 'Flock barycenter', 'Dog', 'Location', 'northwest', ...
    'FontName', font_name, 'FontSize', font_size)
legend('boxoff')

set(gca, 'FontName', font_name, 'FontSize', font_size, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')

axis('equal')
axis([min(min(pos_x))-2 55 min(min(pos_y))-2 max(max(pos_y))+2])

xlabel('x (m)', 'FontName', font_name, 'FontSize', font_size)
ylabel('y (m)', 'FontName', font_name, 'FontSize', font_size)
