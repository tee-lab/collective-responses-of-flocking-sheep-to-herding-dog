%% Run this code to reproduce Fig 7e of the main text

close all
clear
clc

%% Load data

load('hm_n_14.mat')

tm_delay = 20; % max delay till which lags are calculated (it will be 20*0.1 = 2 sec)
cmin = 0.5; % minimum cross-correlation coefficient for it to be considered significant.
no_shp_dg = no_shp + 1;
no_shp_dg = no_shp;

%% Constructing leader-follower network for each drive. 

dij = []; % (xj - xi) \cdot vgroup

node_indeg_final = nan(no_shp_dg, no_it);
di_final = nan(no_shp_dg, no_it);

for i = 1:no_it

    v_s = vel_s(:,:,:,i); % sheep vel
    v_dog = vel_d(:,:,i); % dog vel
    v_d = zeros(1,2,n_iter);
    v_d(1,1,:) = v_dog(:,1);
    v_d(1,2,:) = v_dog(:,2);

    % concatenating sheep and dog velocities to handle the data well.
    vel_shp_dog = cat(1, v_s, v_d);
    
    % vx and vy
    vx = squeeze(vel_shp_dog(:,1,:));
    vy = squeeze(vel_shp_dog(:,2,:));
    vx = vx(:,2:end);
    vy = vy(:,2:end);

    % grp velocity
    vx_grp = mean(vx(1:no_shp,:),1); % grp speed of sheep, x
    vy_grp = mean(vy(1:no_shp,:),1); % grp speed of sheep, y
    
    dr_t = length(vx_grp);
    vel_grp = zeros(2,dr_t);
    vel_grp(1,:) = vx_grp;
    vel_grp(2,:) = vy_grp;
    vel_grp = vel_grp./vecnorm(vel_grp,2,1); 
    
    % concatenating positions of sheep and dog
    pos_shp = pos_s(:,:,:,i);
    pos_dog = pos_d(:,:,i);
    x_d = zeros(1,2,n_iter);
    x_d(1,1,:) = pos_dog(:,1);
    x_d(1,2,:) = pos_dog(:,2);
    pos_temp = cat(1, pos_shp, x_d);
    pos_temp = pos_temp(:,:,2:end);
    
    % storing relative positions
    dij_temp = nan(no_shp_dg,no_shp_dg,dr_t);

    agent_i = []; % leader
    agent_j = []; % follower
    wghts = []; % lag is considered as weights (Nagy et al)

    ccf_lag_temp = nan(no_shp_dg, no_shp_dg); % lag between each pair

    for ind = 1:no_shp_dg

        vx_i = vx(ind,:);
        vy_i = vy(ind,:);

        for j = 1:no_shp_dg

            if j ~= ind

                dijt = squeeze(pos_temp(ind,:,:)) - squeeze(pos_temp(j,:,:));
                dijt = dot(dijt, vel_grp, 1);
                dij_temp(ind,j,:) = dijt;

                vx_j = vx(j,:);
                vy_j = vy(j,:);

                % calculating velocity cross-correlations as defined by
                % Nagy et al
                [ccf_vix_vjx_temp, ~] = xcorr(vx_i, vx_j, tm_delay, 'unbiased'); % x component
                [ccf_viy_vjy_temp, tlag] = xcorr(vy_i, vy_j, tm_delay, 'unbiased'); % y component

                ccf_vel_temp = (ccf_vix_vjx_temp + ccf_viy_vjy_temp); % x + y
                ccf_vel_temp_abs = abs(ccf_vel_temp); % magnitude of ccf
                tlag_id = find(ccf_vel_temp_abs == max(ccf_vel_temp_abs)); % time lag at max ccf

                % check if max(abs(ccf)) is greater than cmin and that the time lag is < 0
                % then i will be leader
                if max(ccf_vel_temp_abs) >= cmin && tlag(tlag_id) < 0

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

        struct_graph = digraph(agent_j, agent_i, wghts, no_shp_dg);
        node_indeg_temp = indegree(struct_graph)';

        mean_dij = mean(dij_temp, 3); % mean distance between i and j
        mean_di = mean(mean_dij, 2, 'omitnan');

        node_indeg_final(:,i) = node_indeg_temp;
        di_final(:,i) = mean_di;

    end

end

%% Plotting mean d_i as a funtion of indegree

% indegree of an individual and its corresponding d_i
node_indeg_final = node_indeg_final(:);
di_final = di_final(:);

max_node_indeg_final = max(node_indeg_final);
di_mean_final = nan(1,max_node_indeg_final+1); % mean di for a given indegree
di_std_err_final = nan(1,max_node_indeg_final+1); % std err for a given indegree
sum_id = nan(1,max_node_indeg_final+1);

for i = 0:max_node_indeg_final
    
    id_sd = find(node_indeg_final == i); % identify indices of a i indegree
    sum_id(i+1) = length(id_sd);
    di_mean_final(1,i+1) = mean(di_final(id_sd)); % calculate mean d_i
    di_std_err_final(1,i+1) = std(di_final(id_sd))/sqrt(length(id_sd));

end

indeg_mean_final = 0:1:max_node_indeg_final;

lin_fit_indig_di_all = polyfit(indeg_mean_final, di_mean_final, 1);
fit_di_all = polyval(lin_fit_indig_di_all, indeg_mean_final);
[rho_all, pval_all] = corr(indeg_mean_final', di_mean_final', 'type', 'Spearman');

fprintf('Pearson correlation for indegree versus d_{i} in model is %.2f and P = %.10f\n', rho_all, pval_all)

fs = 20;
fname = 'Arial';
model_fig_indeg_di = figure('Position', [300 300 800 600]); 
errorbar(indeg_mean_final, di_mean_final, di_std_err_final, 'o', 'color', "#3182bd",...
    'linewidth', 2, 'MarkerSize',15,'MarkerEdgeColor', "#3182bd",'MarkerFaceColor', "#3182bd")
set(gca, 'XLim', [-0.35 13.25], 'XTick', 0:2:13, 'YLim', [-0.5 1.05], 'YTick', -0.4:0.4:1, ...
    'FontName', fname, 'FontSize', fs, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xlabel('Indegree', 'FontName', fname, 'FontSize', fs)
ylabel('$\left<d_{i}\right>$ (m)', 'Interpreter','latex', ...
    'FontName', fname, 'FontSize', fs)
exportgraphics(model_fig_indeg_di, 'model_fig_indeg_di.pdf', 'ContentType', 'vector')