%% Run this code to reproduce Supplementary Fig 8 and 11

close all
clear
clc

%% Loading data

load('drives_data.mat')
load('sheep_all_dat.mat') % loading all sheep data --- pos, vel, etc
tm_delay = 20; % max delay till which lags are calculated (it will be 20*0.1 = 2 sec)
no_shp_dg = no_ind - 2; % total number of individuals considered in this analysis (sheep + dog)
no_shp = no_shp_dg - 1; % no.of sheep
% no_shp_dg = no_shp; % comment this if you want to include dog in the analysis
cmin = 0.5; % minimum cross-correlation coefficient for it to be considered significant.

%% Constructing leader-follower dynamics when speed is considered

node_indeg = []; % storing indegrees
ccf_lag = []; % storing lag time

for ev = 1:length(events)

    evt = events(ev);
    ev_st = eval(strcat('ev_st_', num2str(evt))); % start time 
    ev_et = eval(strcat('ev_et_', num2str(evt))); % end time
    drvs = length(ev_st); % no.of drives in the event

    vel = eval(strcat('vel_ev_',num2str(evt))); % load heading angles (phi)
    
    node_indeg_temp = zeros(drvs, no_shp_dg);

    for dr = 1:drvs

        ev_st_dr = ev_st(dr); % drive start time
        ev_et_dr = ev_et(dr); % drive end time
        
        vel_sheep = vel(1:no_shp_dg,:,ev_st_dr:ev_et_dr); % velocity 
        spd_sheep = squeeze(sqrt(vel_sheep(:,1,:).^2 + vel_sheep(:,2,:).^2)); % speed
            
        agent_i = []; % storing leaders 
        agent_j = []; % storing followers
        wghts = [];
        
        ccf_lag_temp = nan(no_shp_dg, no_shp_dg); % storing lag time

        for ind = 1:no_shp_dg
            
            si = spd_sheep(ind,:) - mean(spd_sheep(ind,:)); % speed of agent i

            for j = 1:no_shp_dg

                if j ~= ind
                    
                    sj = spd_sheep(j,:) - mean(spd_sheep(j,:)); % speed of agent j
                    
                    % calculating cross-correlation
                    [ccf_si_sj, tlag] = xcorr(si, sj, tm_delay, 'normalized');

                    ccf_si_sj_abs = abs(ccf_si_sj);
                    tlag_id = find(ccf_si_sj_abs == max(ccf_si_sj_abs));

                    if max(ccf_si_sj_abs) >= cmin && tlag(tlag_id) < 0
                        
                        ccf_lag_temp(ind,j) = tlag(tlag_id);
                        ccf_lag_temp(j,ind) = -tlag(tlag_id);
                        agent_i = [agent_i ind];
                        agent_j = [agent_j j];
                        wghts = [wghts abs(tlag(tlag_id))];

                    end

                end

            end

        end
        
        if isempty(agent_i) == 0
            struct_graph = digraph(agent_j, agent_i, wghts, no_shp_dg);
            node_indeg_temp(dr,:) = indegree(struct_graph)';
        end
       
        ccf_lag = cat(3, ccf_lag, ccf_lag_temp);

    end

    node_indeg = [node_indeg; node_indeg_temp];

end

%% Calculate average spatial position along front-back axis

dij = []; % (xi - xj) \cdot vgroup
di_mean_all_drive = []; % saving dij for each drive

for ev = 1:length(events)

    evt = events(ev);
    ev_st = eval(strcat('ev_st_', num2str(evt))); % start time 
    ev_et = eval(strcat('ev_et_', num2str(evt))); % end time
    drvs = length(ev_st); % no.of drives in the event

    phi_temp = eval(strcat('phi_ev_',num2str(evt))); % load heading angles (phi)
    pos = eval(strcat('pos_ev_', num2str(evt))); % load position data

    for dr = 1:drvs

        ev_st_dr = ev_st(dr); % drive start time
        ev_et_dr = ev_et(dr); % drive end time
        
        vx = cos(phi_temp(1:no_shp_dg,ev_st_dr:ev_et_dr)); % x component of the velocity
        vy = sin(phi_temp(1:no_shp_dg,ev_st_dr:ev_et_dr)); % y component of the velocity
        vx = vx(:,2:end); % taking from 2nd point because vx(:,1) is not possible to calculate (backward difference)
        vy = vy(:,2:end); % similar to vx

        vx_grp = mean(vx(1:no_shp,:),1); % grp speed of sheep, x
        vy_grp = mean(vy(1:no_shp,:),1); % grp speed of sheep, y
        dr_t = length(vx_grp); % length of the data
        vel_grp = zeros(2,dr_t); % group velocity
        vel_grp(1,:) = vx_grp;
        vel_grp(2,:) = vy_grp;
        vel_grp = vel_grp./vecnorm(vel_grp,2,1); 
        pos_temp = pos(1:no_shp_dg,:,ev_st_dr:ev_et_dr); % position of sheep + dog
        pos_temp = pos_temp(:,:,2:end);

        dij_temp = nan(no_shp_dg,no_shp_dg,dr_t); % dij for a given drive

        for ind = 1:no_shp_dg

            for j = 1:no_shp_dg

                if j ~= ind

                    dijt = squeeze(pos_temp(ind,:,:)) - squeeze(pos_temp(j,:,:)); % xi - xj (dimension is 2 x t)
                    dijt = dot(dijt, vel_grp, 1); % (xi - xj) \cdot vg
                    dij_temp(ind,j,:) = dijt;

                end

            end

        end
        
        dij = cat(3, dij, dij_temp); % concatanate them.
        mean_dij_temp = mean(dij_temp, 3); % mean distance between i and j
        mean_di_temp = mean(mean_dij_temp, 2, 'omitnan'); % mean position of i
        di_mean_all_drive = [di_mean_all_drive; mean_di_temp'];

    end

end

mean_dij = mean(dij, 3); % mean distance between i and j across all drives
mean_di = mean(mean_dij, 2, 'omitnan'); % di (avg over j)
std_di = std(mean_dij, 0, 2, 'omitnan')/sqrt(no_shp_dg);

%% Figure for mean dij vs sheep (Supplementary Fig 8)

figure_dij_shp = figure('Position', [300 300 900 700]);
font_name = 'Arial';

[~, ord] = sort(mean_di, 'descend');

for i = 1:no_shp_dg
    
    ord_temp = ord(i);
    plot(i*ones(no_shp_dg,1), mean_dij(ord_temp,:), 'o', 'Color', '#C3C3C3', ...
        'MarkerSize', 5, 'MarkerFaceColor', '#C3C3C3')
    hold on

end

errorbar(1:no_shp_dg, mean_di(ord), std_di, 'o', 'color', 'k',...
    'linewidth', 2, 'MarkerSize',15,'MarkerEdgeColor', "k", ...
    'MarkerFaceColor', "k")
set(gca, 'XLim', [0.5 14.2], 'XTick', 1:14, 'XTickLabel', string(ord),...
    'YLim', [-1.5 1.5], 'YTick', -1.5:0.5:1.5, 'FontName', font_name, ...
    'FontSize', 30, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Sheep ID', 'FontSize', 30, 'FontName', ...
    font_name, 'Color', 'k')
ylabel('$d_{ij}$ (m)', 'Interpreter','latex', ...
    'FontName', font_name, 'FontSize', 30, 'Color', 'k')
% exportgraphics(figure_dij_shp, 'dij_sheep.pdf', 'ContentType', 'vector')

%% indeg vs dij similar to analysis done in herding model (Supplementary Fig 11)

di_mean_all_drive = di_mean_all_drive';
node_indeg = node_indeg';
di_mean_all_drive = di_mean_all_drive(:);
node_indeg = node_indeg(:);

di_all_drive_mean = nan(1,max(node_indeg)+1); % all possible indegree
di_all_drive_std_err = nan(1,max(node_indeg)+1); % mean di corresponding to given indegree
sum_id = nan(1,max(node_indeg)+1);

for i = 0:10%max(node_indeg)

    id_sd = find(node_indeg == i);
    sum_id(i+1) = length(id_sd); % no.of times we observe given indegree
    di_all_drive_mean(1,i+1) = mean(di_mean_all_drive(id_sd));
    di_all_drive_std_err(1,i+1) = std(di_mean_all_drive(id_sd))/sqrt(length(id_sd));

end

% indeg_all_drive_mean = 0:1:max(node_indeg);
indeg_all_drive_mean = 0:1:10;
di_all_drive_mean = di_all_drive_mean(~isnan(di_all_drive_mean));
di_all_drive_std_err = di_all_drive_std_err(~isnan(di_all_drive_std_err));

[rho_all, pval_all] = corr(indeg_all_drive_mean', di_all_drive_mean', 'type', 'Pearson');
fprintf('Pearson correlation for indegree versus d_{i} in sheep herd is %.2f and P = %.2f\n', rho_all, pval_all)

fig_dij_indeg = figure('Position', [300 300 900 700]);

plot(indeg_all_drive_mean, di_all_drive_mean, '^', 'color', "#3182bd",...
    'MarkerSize',18,'MarkerEdgeColor', "#3182bd",'MarkerFaceColor', "#3182bd")
set(gca, 'XLim', [-0.2 10.2], 'XTick', 0:2:13, 'YLim', [-1.2 0.8] ,'YTick', -1.2:0.4:0.8, ...
    'FontName', font_name, 'FontSize', 30, 'XColor', 'k', 'YColor', 'k')
xlabel('Indegree','FontName', font_name, 'FontSize', 30, 'Color', 'k')
ylabel('$\left<d_{i}\right>$ (m)', 'Interpreter','latex', ...
    'FontName', font_name, 'FontSize', 30, 'Color', 'k')
% exportgraphics(fig_dij_indeg, 'sm_net_speed.pdf', 'ContentType', 'vector')