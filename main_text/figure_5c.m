close all
clear
clc

%% Calculation for drives only
load('drives_data.mat') % can load both 'ext_drives_data.mat' and 'drives_data.mat'
load('sheep_all_dat.mat') % loading all sheep data --- pos, vel, etc
tm_delay = 20; % max delay till which lags are calculated (it will be 20*0.1 = 2 sec)
no_shp_dg = no_ind - 1; % total number of individuals considered in this analysis (sheep + dog)
no_shp = no_shp_dg - 1; % no.of sheep
% no_shp_dg = no_shp;
cmin = 0.5; % minimum cross-correlation coefficient for it to be considered significant. 
fs = 40; % font size

%% Constructing leader-follower network for each drive. 

%%%% NOTE: To obtain Figure 5c1 and 4c2 set ev = 2 and dr = 3;

fig_net = 1; % network for a given drive
fig_ccf = 1; % cc function for a given drive
ccf_lag = []; % lag for given pair of individuals
dr_t = 0;
no_dir_edges = 0;
tot_drvs = numel(ev_et_1) + numel(ev_et_2) + numel(ev_et_3) + numel(ev_et_4) + ...
    numel(ev_et_5) + numel(ev_et_6) + numel(ev_et_8) + numel(ev_et_9);
tot_no_dir_edges = tot_drvs*(factorial(no_shp_dg)/(factorial(no_shp_dg - 2)*factorial(2)));

node_label = 1:no_shp_dg; % node label - either sheep (1-14) or dog
dg_index = node_label == 15; % dog index
node_label = string(node_label); % making the sheep index into string
node_label(dg_index) = "Dog";

for ev = 1:length(events)

    evt = events(ev); % event
    ev_st = eval(strcat('ev_st_', num2str(evt))); % start time of all drives in a given event
    ev_et = eval(strcat('ev_et_', num2str(evt))); % end time of all drives in a given event
    drvs = length(ev_st); % no.of drives in the event

    phi_temp = eval(strcat('phi_ev_',num2str(evt))); % load heading angles (phi)
    pos = eval(strcat('pos_ev_',num2str(evt))); % load position
    vel = eval(strcat('vel_ev_',num2str(evt))); % load velocity

    node_indeg_temp = zeros(drvs, no_shp_dg); % no-of indegrees of individual in each drive

    for dr = 1:drvs
        
        dr_t = dr_t + 1;
        ev_st_dr = ev_st(dr); % drive start time
        ev_et_dr = ev_et(dr); % drive end time
        
        vx = cos(phi_temp(1:no_shp_dg,ev_st_dr:ev_et_dr)); % x coor of vel
        vy = sin(phi_temp(1:no_shp_dg,ev_st_dr:ev_et_dr)); % y coor of vel
        vx = vx(:,2:end);
        vy = vy(:,2:end);

        pos_sheep = pos(1:no_shp,:,ev_st_dr:ev_et_dr); % sheep position
        vel_sheep = vel(1:no_shp,:,ev_st_dr:ev_et_dr); % velocity of sheep
        pos_shp_dg = pos(:,:,ev_st_dr:ev_et_dr);

        grp_vel = squeeze(mean(vel_sheep, 1)); % grp velocity
        grp_pos = mean(pos_sheep, 1); % grp centre

        pos_sheep_dog_cfr = pos_shp_dg(1:no_shp_dg,:,:) - grp_pos; % position relative to centre

        tm = size(vel_sheep,3); % length of the data set
        rel_x_temp = nan(no_shp_dg,tm); % relative x position
        rel_y_temp = nan(no_shp_dg,tm); % relative y position

        for t = 1:tm

            rot_ang = atan2(grp_vel(2,t), grp_vel(1,t)); % angle made by grp velocity along x-axis
            rot_mat = [cos(rot_ang) sin(rot_ang); -sin(rot_ang) cos(rot_ang)]; % rotation matrix towards x axis
            pos_shp_cfr_temp = pos_sheep_dog_cfr(:,:,t)'; % the vector should be of the for dimension X no.of particles
            new_pos = rot_mat*pos_shp_cfr_temp; % rotated coordinates
            new_pos = [0 -1; 1 0]*new_pos; % rotate towards y axis
            new_pos = new_pos'; % new positon
            rel_x_temp(:,t) = new_pos(:,1);
            rel_y_temp(:,t) = new_pos(:,2);

        end

        rel_x = mean(rel_x_temp,2);
        rel_y = mean(rel_y_temp,2);

        agent_i = []; % leader 
        agent_j = []; % follower
        wghts = []; % lag is considered as weights (Nagy et al)
        
        ccf_lag_temp = nan(no_shp_dg, no_shp_dg); % lag between each pair

        for ind = 1:no_shp_dg

            vx_i = vx(ind,:); % vel_x of i
            vy_i = vy(ind,:); % vel_y of i

            for j = 1:no_shp_dg

                if j ~= ind

                    vx_j = vx(j,:); % vel_x of j
                    vy_j = vy(j,:); % vel_y of j
                    
                    % calculating velocity cross-correlations as defined by
                    % Nagy et al
                    [ccf_vix_vjx_temp, ~] = xcorr(vx_i, vx_j, tm_delay, 'unbiased'); % x component
                    [ccf_viy_vjy_temp, tlag] = xcorr(vy_i, vy_j, tm_delay, 'unbiased'); % y component

                    ccf_vel_temp = (ccf_vix_vjx_temp + ccf_viy_vjy_temp); % x + y
                    ccf_vel_temp_abs = abs(ccf_vel_temp); % magnitude of ccf
                    % tlag_id = find(ccf_vel_temp_abs == max(ccf_vel_temp_abs)); % time lag at max ccf
                    tlag_id = find(ccf_vel_temp == max(ccf_vel_temp));

                    % if max(ccf_vel_temp_abs) >= cmin && tlag(tlag_id) < 0
                    if max(ccf_vel_temp) >= cmin && tlag(tlag_id) < 0
                        
                        no_dir_edges = no_dir_edges + 1;
                        ccf_lag_temp(ind,j) = tlag(tlag_id);
                        ccf_lag_temp(j,ind) = -tlag(tlag_id);
                        agent_i = [agent_i ind];
                        agent_j = [agent_j j];
                        wghts = [wghts abs(tlag(tlag_id))];

                    end

                end

            end

        end
        
        % if there are some leader-follower then, calculate indegree for
        % each of them
        
        if isempty(agent_i) == 0
            
            fs = 30;
            struct_graph = digraph(agent_j, agent_i, wghts, no_shp_dg); % construct the directed network
            node_indeg_temp(dr,:) = indegree(struct_graph)';

            fig_net_dog = figure('Position', [300 300 800 600]);
            plot(struct_graph, 'XData', rel_x, 'YData', rel_y, 'Marker', 'o', 'NodeColor', '#2c7fb8', ...
                'MarkerSize', 1.5*(indegree(struct_graph)'+1), 'LineWidth', 2, 'NodeLabel', node_label,...
                'EdgeColor', '#bcbddc', 'ArrowSize', 8, 'NodeFontSize', fs, 'NodeFontName', 'Arial')
            hold on
            q = quiver(0, 0, 0, 1, 0.6, 'LineWidth', 3);
            q.Marker = 'none'; q.Color = '#d95f0e';
            hold on
            plot(0, 0, '.', 'Color', '#d95f0e', 'MarkerSize', 40)

            set(gca, 'YTick', -5:2:1, 'FontName', 'Arial', 'FontSize', fs)
            % set(gca, 'YTick', -0.9:.3:1, 'FontName', 'Arial', 'FontSize', fs)
            xlabel('x (m)', 'FontName', 'Arial', 'FontSize', fs)
            ylabel('y (m)', 'FontName', 'Arial', 'FontSize', fs)

            % title(strcat('Video- ', num2str(evt), ' and Drive-  ', num2str(dr)))

        end

        fig_net = fig_net + 1;
        ccf_lag = cat(3, ccf_lag, ccf_lag_temp); % concatenating ccf lags

    end
    
    fig_ccf = fig_ccf + 1;

end