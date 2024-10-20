%% Run this code to reproduce Fig.7d-e of the main text. 
    
close all
clear
clc
    
%% Calculation for all drives
load('ext_drives_data.mat') % 
load('sheep_all_dat.mat') % loading all sheep data --- pos, vel, etc
tm_delay = 20; % max delay till which lags are calculated (it will be 20*0.1 = 2 sec)
no_shp_dg = no_ind - 1; % total number of individuals considered in this analysis (sheep + dog)
no_shp = no_shp_dg - 1; % no.of sheep
no_shp_dg = no_shp; % comment this code when you want to add dog to the analysis
cmin = 0.5; % minimum cross-correlation coefficient for it to be considered significant. 

font_name = 'Arial';
font_size = 25;

%% Constructing leader-follower network for each drive. 

node_indeg = []; % node indegree (number of followers a given node has)
ccf_lag = []; % lag for given pair of individuals

for ev = 1:length(events)

    evt = events(ev); % event
    ev_st = eval(strcat('ev_st_', num2str(evt))); % start time of all drives in a given event
    ev_et = eval(strcat('ev_et_', num2str(evt))); % end time of all drives in a given event
    drvs = length(ev_st); % no.of drives in the event

    phi_temp = eval(strcat('phi_ev_',num2str(evt))); % load heading angles (phi)
    node_indeg_temp = zeros(drvs, no_shp_dg); % no-of indegrees of individual in each drive

    for dr = 1:drvs
        
        ev_st_dr = ev_st(dr); % drive start time
        ev_et_dr = ev_et(dr); % drive end time
        
        vx = cos(phi_temp(1:no_shp_dg,ev_st_dr:ev_et_dr)); % vel_x for a given drive
        vy = sin(phi_temp(1:no_shp_dg,ev_st_dr:ev_et_dr)); % vel_y for a given drive
        vx = vx(:,2:end); 
        vy = vy(:,2:end);
            
        agent_i = []; % store index of all leaders
        agent_j = []; % store index of all followers
        wghts = []; % lag is considered as weights (Nagy et al)
        
        ccf_lag_temp = nan(no_shp_dg, no_shp_dg); % lag between each pair
        
        % calculate cross-correlation for all pairs
        for ind = 1:no_shp_dg

            vx_i = vx(ind,:);
            vy_i = vy(ind,:);

            for j = 1:no_shp_dg

                if j ~= ind

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
            node_indeg_temp(dr,:) = indegree(struct_graph)';
            
        end
       
        ccf_lag = cat(3, ccf_lag, ccf_lag_temp); % concatenating ccf lags

    end

    node_indeg = [node_indeg; node_indeg_temp];

end

node_indeg_mean = mean(node_indeg,1); % mean of indegrees over all drives and events
[~, ord] = sort(node_indeg_mean, 'descend'); % ordering from agent with highest indegree to the lowest

%% Calculate relative distance along front-back axis

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
                    dijt = dot(dijt, vel_grp, 1); % (xi - xj) \cdot v_B
                    dij_temp(ind,j,:) = dijt;

                end

            end

        end
        
        dij = cat(3, dij, dij_temp); % concatanate them.
        mean_dij_temp = mean(dij_temp, 3); % mean distance between i and j
        mean_di_temp = mean(mean_dij_temp, 2, 'omitnan'); % mean position of i in the flock
        di_mean_all_drive = [di_mean_all_drive; mean_di_temp'];

    end

end

mean_dij = mean(dij, 3); % mean distance between i and j
mean_di = mean(mean_dij, 2, 'omitnan'); 
std_di = std(mean_dij, 0, 2, 'omitnan')/sqrt(no_shp_dg);

%% indeg vs dij 

di_mean_all_drive = di_mean_all_drive';
node_indeg = node_indeg';
di_mean_all_drive = di_mean_all_drive(:);
node_indeg = node_indeg(:);

di_all_drive_mean = nan(1,max(node_indeg)+1);
di_all_drive_std_err = nan(1,max(node_indeg)+1);
sum_id = nan(1,max(node_indeg)+1);
exp_less_zero = nan(1,max(node_indeg)+1);
exp_more_zero = nan(1,max(node_indeg)+1);

for i = 0:max(node_indeg)

    id_sd = find(node_indeg == i);
    sum_id(i+1) = length(id_sd);
    di_all_drive_mean(1,i+1) = mean(di_mean_all_drive(id_sd));
    di_all_drive_std_err(1,i+1) = std(di_mean_all_drive(id_sd))/sqrt(length(id_sd));
    exp_less_zero(i+1) = sum(di_mean_all_drive(id_sd) <= 0)/length(di_mean_all_drive(id_sd));
    exp_more_zero(i+1) = sum(di_mean_all_drive(id_sd) > 0)/length(di_mean_all_drive(id_sd));

end

indeg_all_drive_mean = 0:1:max(node_indeg);

% lin_fit_indig_di_all_drive = polyfit(indeg_all_drive_mean, di_all_drive_mean, 1);
% fit_di_all_drive = polyval(lin_fit_indig_di_all_drive, indeg_all_drive_mean);
[rho_all, pval_all] = corr(indeg_all_drive_mean', di_all_drive_mean', 'type', 'Pearson');
fprintf('Pearson correlation for indegree versus d_{i} in sheep herd is %.2f and P = %.10f\n', rho_all, pval_all)

%% Plotting Fig.7d of main text

% fig_7d = figure('Position', [300 300 800 600]);
fig_all_data = figure('Position', [300 300 1400 450]);

node_indeg_str = [];
di_mean_bp = [];

subplot(1,2,1)
for i = 0:max(node_indeg)

    id_sd = find(node_indeg == i);
    di_indeg = di_mean_all_drive(id_sd);
    di_mean_bp = [di_mean_bp; di_indeg];
    node_indeg_str = [node_indeg_str; i*ones(length(id_sd),1)];
    plot(i*ones(length(id_sd),1), di_indeg, 'o', 'Color', "#93C0DF", 'MarkerSize', 5, ...
        'MarkerFaceColor', '#93C0DF')
    hold on
    
end

% ylim([-2.8 2.8])
 
% node_indeg_str = num2str(node_indeg_str);
% boxplot(di_mean_bp, node_indeg_str, 'OutlierSize', 5, 'Notch', 'on', 'Whisker', 1)
% ylim([-2.04 1.78])

errorbar(indeg_all_drive_mean, di_all_drive_mean, di_all_drive_std_err, 'o', 'color', "#3182bd",...
    'linewidth', 2, 'MarkerSize',10,'MarkerEdgeColor', "#3182bd",'MarkerFaceColor', "#3182bd")
% plot(indeg_all_drive_mean, di_all_drive_mean, '^', 'color', "#3182bd",...
%     'MarkerSize',18,'MarkerEdgeColor', "#3182bd",'MarkerFaceColor', "#3182bd")
set(gca, 'XLim', [-0.35 13.25], 'XTick', 0:2:13, 'YLim', [min(di_mean_all_drive) max(di_mean_all_drive)], ...
    'YTick', -4:1:5, 'LineWidth', 1, 'FontName', font_name, 'FontSize', ...
    font_size, 'Xcolor', 'k', 'YColor', 'k')
% set(gca, 'XLim', [-0.35 13.25], 'XTick', 0:2:13, 'YLim', [-1 0.9], ...
%     'YTick', -0.9:.3:1, 'LineWidth', 1, 'FontName', font_name, 'FontSize', ...
%     font_size, 'Xcolor', 'k', 'YColor', 'k')
xlabel('Indegree', 'FontName', font_name, 'FontSize', font_size);
ylabel('$d_{i}$ (m)', 'Interpreter','latex', ...
    'FontName', font_name, 'FontSize', 40);
% ylabel('$\left<d_{i}\right>$ (m)', 'Interpreter','latex', 'FontName', font_name, 'FontSize', 40);
% exportgraphics(fig_7d, 'figure_7d_all_data.pdf', 'ContentType', 'vector')

%% Plotting Fig.7e of main text

% Data is give in the file, 'model_di_indeg.csv'

model_di_indeg = readmatrix('model_di_indeg.csv');
indeg_mean_final = model_di_indeg(1,:);
di_mean_final = model_di_indeg(2,:);
di_std_err_final = model_di_indeg(3,:);

% lin_fit_indig_di_all = polyfit(indeg_mean_final, di_mean_final, 1);
% fit_di_all = polyval(lin_fit_indig_di_all, indeg_mean_final);
[rho_all, pval_all] = corr(indeg_mean_final', di_mean_final', 'type', 'Spearman');

fprintf('Pearson correlation for indegree versus d_{i} in model is %.2f and P = %.10f\n', rho_all, pval_all)

model_di_indeg_all_data = readmatrix('model_di_indeg_all_data.csv');
model_node_indeg = model_di_indeg_all_data(:,1);
model_di = model_di_indeg_all_data(:,2);

fig_7e = figure('Position', [300 300 800 600]);

model_node_indeg_str = [];
model_di_mean_bp = [];
% model_di_mean = nan(length(model_di),length(indeg_mean_final));
model_sum_id = nan(1,max(model_node_indeg)+1);
model_less_zero = nan(1,max(model_node_indeg)+1);
model_more_zero = nan(1,max(model_node_indeg)+1);

% subplot(1,2,2)
for i = 0:max(model_node_indeg)

    id_sd = find(model_node_indeg == i);
    di_indeg = model_di(id_sd);
    model_sum_id(i+1) = length(id_sd);
    model_less_zero(i+1) = sum(di_indeg <= 0)/length(di_indeg);
    model_more_zero(i+1) = sum(di_indeg > 0)/length(di_indeg);
    model_di_mean_bp = [model_di_mean_bp; di_indeg];
    % model_di_mean(1:length(id_sd),i+1) = di_indeg;
    model_node_indeg_str = [model_node_indeg_str; i*ones(length(id_sd),1)];
    % plot(i*ones(length(id_sd),1), di_indeg, 'o', 'Color', "#93C0DF", 'MarkerSize', 3, ...
    %     'MarkerFaceColor', '#93C0DF')
    % hold on
    
end
 
% model_node_indeg_str = num2str(model_node_indeg_str);
% boxplot(model_di_mean_bp, model_node_indeg_str, 'OutlierSize', 0.001, 'Notch', 'on', 'Whisker', 1)
% ylim([-2.8 2.8])


% errorbar(indeg_mean_final, di_mean_final, di_std_err_final, 'o', 'color', "#3182bd",...
%     'linewidth', 2, 'MarkerSize',8,'MarkerEdgeColor', "#3182bd",'MarkerFaceColor', "#3182bd")

plot(indeg_mean_final, di_mean_final, '^', 'color', "#3182bd",...
    'MarkerSize',18,'MarkerEdgeColor', "#3182bd",'MarkerFaceColor', "#3182bd")

set(gca, 'XLim', [-0.35 13.25], 'XTick', 0:2:13, 'YLim', [-1 0.9], ...
    'YTick', -0.9:0.3:0.9, 'LineWidth', 1, 'FontName', font_name, 'FontSize', ...
    font_size, 'Xcolor', 'k', 'YColor', 'k')
% set(gca, 'XLim', [-0.35 13.25], 'XTick', 0:2:13, 'YLim', [-5 5], ...
%     'YTick', -10:1:10, 'LineWidth', 1, 'FontName', font_name, 'FontSize', ...
%     font_size, 'Xcolor', 'k', 'YColor', 'k')
xlabel('Indegree', 'FontName', font_name, 'FontSize', font_size);
% ylabel('$d_{i} $', 'Interpreter','latex', 'FontName', font_name, 'FontSize', 40);
ylabel('$\left<d_{i}\right>$', 'Interpreter', 'latex', 'FontName', font_name, 'FontSize', 40);

exportgraphics(fig_7e, 'figure_7e.pdf', 'ContentType', 'vector')
% exportgraphics(fig_all_data, 'figure_7_all_data.pdf', 'ContentType', 'vector')

%% proportion front and back

fig_7f = figure('Position', [300 300 800 600]);

plot(indeg_all_drive_mean, exp_more_zero, 'o', 'color', "#3182bd",'MarkerSize', 15, ...
    'MarkerFaceColor', '#3182bd')
hold on
plot(indeg_mean_final, model_more_zero, 's', 'Color', "#93C0DF",'MarkerSize', 15, ...
    'MarkerFaceColor', '#93C0DF')

set(gca, 'XLim', [-0.35 13.25], 'XTick', 0:2:13, 'YLim', [0 1], ...
    'YTick', 0:.2:1, 'LineWidth', 1, 'FontName', font_name, 'FontSize', ...
    font_size, 'Xcolor', 'k', 'YColor', 'k')
xlabel('Indegree', 'FontName', font_name, 'FontSize', font_size);
% yl = ylabel('$\left<d_{i}\right>$', 'Interpreter','latex', ...
%     'FontName', font_name, 'FontSize', 40);
ylabel('Proportion of time d_i > 0', 'FontName', font_name, 'FontSize', font_size);

legend({'Data', 'Model'}, 'Location', 'best', 'FontName', font_name, 'FontSize', font_size)
legend('boxoff')

[r,p] = corr(indeg_all_drive_mean', exp_more_zero', 'Type', 'Pearson')
[rm, pm] = corr(indeg_mean_final', model_more_zero', 'Type', 'Pearson')

exportgraphics(fig_7f, 'figure_7f.pdf', 'ContentType', 'vector')