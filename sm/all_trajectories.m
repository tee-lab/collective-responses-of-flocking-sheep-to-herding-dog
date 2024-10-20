%% plotting all trajectories to check if they spend a lot of time near boundary. 

close all
clear
clc         

%% Loading data

load('sheep_all_dat.mat') 
lw_plt = 1.5; % width of plot line 
lw_axis = 1.5; % width of axis
fn = 'Arial'; % font name
fs = 20; % font size
fig_pos = [300 300 600 500];
tar_rad = 6; % target radius
fence_col = [0.5 0.5 0.5]; % fence color
lw_fence = 2; % fence line width
lw_tar = 3; % line width of target
no_shp_dg = no_ind - 1; % no.of sheep + dog

%% Plotting trajectories.

fig_all_dat = figure('Position', fig_pos);

for e = 1:no_evnts
    pos_temp = eval(strcat('pos_ev_',num2str(e)));
    for i = 1:no_shp_dg
        if i == 1
            h1 = plot(squeeze(pos_temp(i,1,:)), squeeze(pos_temp(i,2,:)), 'Color', '#56B4E9', ...
            'LineWidth', lw_plt);
            hold on
        elseif i == 15
            h2 = plot(squeeze(pos_temp(i,1,:)), squeeze(pos_temp(i,2,:)), 'Color', 'r', ...
            'LineWidth', lw_plt);
        else
            plot(squeeze(pos_temp(i,1,:)), squeeze(pos_temp(i,2,:)), 'Color', '#56B4E9', ...
            'LineWidth', lw_plt)
        end
    end
end

% plotting fence boundary. coordinates are given below.
plot([0 0], [-25 25], 'Color', fence_col, 'LineWidth', lw_fence)
hold on
plot([0 80], [25 25], 'Color', fence_col, 'LineWidth', lw_fence)
hold on
plot([80 80], [-25 25], 'Color', fence_col, 'LineWidth', lw_fence)
hold on
plot([0 80], [-25 -25], 'Color', fence_col, 'LineWidth', lw_fence)

% targets
h3 = viscircles([10 0], tar_rad, 'Color', '#964B00', 'LineWidth', lw_tar);
hold on
viscircles([70 0], tar_rad, 'Color', '#964B00', 'LineWidth', lw_tar)
hold off

axis('equal')
set(gca, 'XLim', [-10 105], 'XTick', 0:20:80, 'YLim', [-35 35], 'YTick', -25:10:25, ...
    'FontName', fn, 'FontSize', fs, 'LineWidth', lw_axis, 'XColor', 'k', 'YColor', 'k')
xlabel('x (m)', 'FontSize', fs, 'FontName', fn)
ylabel('y (m)', 'FontSize', fs, 'FontName', fn)

lg_plt = [h1 h2 h3];
legend(lg_plt, 'Sheep', 'Dog', 'Target', 'Location', 'east', 'FontSize', fs-8)
legend('boxoff')
exportgraphics(fig_all_dat, 'traj_all_data.jpg', 'Resolution', 600)

%% Plotting trajectories of only drives. 

load('drives_data.mat') % load drives data
load('sheep_all_dat.mat') % loading all sheep data --- pos, vel, etc

%% Reproduce Supplementary Fig 14.

fig_all_drive = figure('Position', fig_pos);

for ev = 1:length(events)

    evt = events(ev); % event
    ev_st = eval(strcat('ev_st_', num2str(evt))); % start time of all drives in a given event
    ev_et = eval(strcat('ev_et_', num2str(evt))); % end time of all drives in a given event
    drvs = length(ev_st); % no.of drives in the event

    pos_temp = eval(strcat('pos_ev_',num2str(evt)));

    for dr = 1:drvs

        ev_st_dr = ev_st(dr); % drive start time
        ev_et_dr = ev_et(dr); % drive end time
        pos_temp_dr = pos_temp(:,:,ev_st_dr:ev_et_dr);

        for i = 1:no_shp_dg
            if i == 1
                h1 = plot(squeeze(pos_temp_dr(i,1,:)), squeeze(pos_temp_dr(i,2,:)), 'Color', '#56B4E9', ...
                    'LineWidth', lw_plt);
                hold on
            elseif i == 15

                h2 = plot(squeeze(pos_temp_dr(i,1,:)), squeeze(pos_temp_dr(i,2,:)), 'Color', 'r', ...
                    'LineWidth', lw_plt);
            else
                plot(squeeze(pos_temp_dr(i,1,:)), squeeze(pos_temp_dr(i,2,:)), 'Color', '#56B4E9', ...
                    'LineWidth', lw_plt)
            end
        end

    end

end

plot([0 0], [-25 25], 'Color', fence_col, 'LineWidth', lw_fence)
hold on
plot([0 80], [25 25], 'Color', fence_col, 'LineWidth', lw_fence)
hold on
plot([80 80], [-25 25], 'Color', fence_col, 'LineWidth', lw_fence)
hold on
plot([0 80], [-25 -25], 'Color', fence_col, 'LineWidth', lw_fence)

% targets
h3 = viscircles([10 0], tar_rad, 'Color', '#964B00', 'LineWidth', lw_tar);
hold on
viscircles([70 0], tar_rad, 'Color', '#964B00', 'LineWidth', lw_tar)
hold off

axis('equal')
set(gca, 'XLim', [-10 105], 'XTick', 0:20:80, 'YLim', [-35 35], 'YTick', -25:10:25, ...
    'FontName', fn, 'FontSize', fs, 'LineWidth', lw_axis, 'XColor', 'k', 'YColor', 'k')
xlabel('x (m)', 'FontSize', fs, 'FontName', fn)
ylabel('y (m)', 'FontSize', fs, 'FontName', fn)

lg_plt = [h1 h2 h3];
legend(lg_plt, 'Sheep', 'Dog', 'Target', 'Location', 'east', 'FontSize', fs-6)
legend('boxoff')
exportgraphics(fig_all_drive, 'traj_all_drive.pdf', 'ContentType', 'vector')