%% Run this code to reproduce Fig.5a-b of the main text. 

close all
clear
clc

%% Loading data

load('ext_drives_data.mat') 
no_evnts = length(events);

no_shp = no_ind - 2; % no.of sheep, i.e., no.of individuals - (dog + shepherd)
font_name = 'Arial';
font_size = 20;

%% Calculate turning angles as a function of viewing angle

del_grp_phi_final = []; % change in group orientation (\Delta \phi_B)
r_bary_dg_ci_final = []; % \psi_BD angle
del_dog_phi_final = []; % \Delta \phi_D of dog
r_dg_bary_ci_final = []; % \psi_{DB} 
del_t = 30; % \Delta t over which turning angle is calculated.

for ev = 1:no_evnts

    pos = eval(strcat('pos_ev_', num2str(ev))); % load position data
    vel = eval(strcat('vel_ev_',num2str(ev))); % load velocity

    vel_sheep = vel(1:14,:,2:end); % velocity of sheep
    grp_vel = squeeze(mean(vel_sheep, 1)); % group velocity
    grp_vel = grp_vel./vecnorm(grp_vel,2,1);
    grp_phi = atan2(grp_vel(2,:), grp_vel(1,:)); % group orientation

    % del_grp_phi = dot(grp_vel(:,1:(end-del_t)), grp_vel(:,del_t+1:end), 1);
    % del_grp_phi = acos(del_grp_phi);
    
    del_grp_phi = nan(length(grp_phi)-del_t,1); % calculate \del\phi = \phi(t+dt) - \phi(t)

    for t = 1:(length(del_grp_phi))

        del_grp_phi_temp = grp_phi(t+del_t) - grp_phi(t); % \Delta \phi of group

        if del_grp_phi_temp > pi 
            del_grp_phi_temp = del_grp_phi_temp - 2*pi;
        elseif del_grp_phi_temp < -pi 
            del_grp_phi_temp = del_grp_phi_temp + 2*pi;
        end

        del_grp_phi(t) = del_grp_phi_temp;

    end
   
    del_grp_phi_final = [del_grp_phi_final del_grp_phi'];

    dog_vel = squeeze(vel(15,:,2:end));
    dog_phi = atan2(dog_vel(2,:), dog_vel(1,:));

    % calculating turning angle of dog
    del_dog_phi = nan(length(dog_phi)-del_t,1);

    for t = 1:length(del_dog_phi)

        del_dog_phi_temp = dog_phi(t+del_t) - dog_phi(t);

        if del_dog_phi_temp > pi 
            del_dog_phi_temp = del_dog_phi_temp - 2*pi;
        elseif del_dog_phi_temp < -pi 
            del_dog_phi_temp = del_dog_phi_temp + 2*pi;
        end

        del_dog_phi(t) = del_dog_phi_temp;
 
    end

    del_dog_phi_final = [del_dog_phi_final del_dog_phi'];

    pos_shp = pos(1:14,:,2:(end-del_t)); % position of sheep 
    vel_shp = vel(1:14,:,2:(end-del_t)); % vel of sheep
    pos_dog = squeeze(pos(15,:,2:(end-del_t))); % pos of dog
    tm = size(pos_shp,3);
    
    r_bary_dg_ci = nan(1,tm); % where is dog from barycentre
    r_dg_bary_ci = nan(1,tm); % where is sheep from dogs prespective
    

    for t = 1:tm

        r_bary_dg_temp = pos_dog(:,t).' - mean(pos_shp(:,:,t),1); % displacement vec between sheep bary and dog
        r_bary_dg_temp = r_bary_dg_temp./sqrt(r_bary_dg_temp(:,1).^2 + r_bary_dg_temp(:,2).^2); % normalising 
        r_dg_bary_temp = -r_bary_dg_temp; % displacement vec between dog and sheep bary

        vel_bary_temp = mean(vel_shp(:,:,t),1); % vel of bary
        vel_bary_temp = vel_bary_temp./vecnorm(vel_bary_temp,2,2); % normalised vel of bary
        vel_dog_temp = dog_vel(:,t); % vel of dog
        vel_dog_temp = vel_dog_temp';
        vel_dog_temp = vel_dog_temp/(vecnorm(vel_dog_temp,2,2)+eps);

        ci_bary_dg_temp = dot(r_bary_dg_temp, vel_bary_temp, 2); 
        ci_bary_dg_temp = acos(ci_bary_dg_temp); % angle btw bary and dog orientation
        ci_bary_dg_sign = cross([vel_bary_temp 0], [r_bary_dg_temp 0], 2);
        ci_bary_dg_temp = ci_bary_dg_temp*sign(ci_bary_dg_sign(1,3)); % sign of the angle
        r_bary_dg_ci(1,t) = ci_bary_dg_temp;  % where is dog from barycentre

        ci_dg_bary_temp = dot(r_dg_bary_temp, vel_dog_temp, 2); % angle btw dog and sheep from dog's orientation
        ci_dg_bary_temp = acos(ci_dg_bary_temp);
        ci_dg_bary_sign = cross([vel_dog_temp 0], [r_dg_bary_temp 0], 2);
        ci_dg_bary_temp = ci_dg_bary_temp*sign(ci_dg_bary_sign(1,3)); % sign og angle
        r_dg_bary_ci(1,t) = ci_dg_bary_temp;

    end

    r_bary_dg_ci_final = [r_bary_dg_ci_final r_bary_dg_ci];
    r_dg_bary_ci_final = [r_dg_bary_ci_final r_dg_bary_ci];

end

%% Turning angle of barycenter 

psi_edges = -pi-pi/12:pi/6:pi+pi/12; % \psi_B edges

[hs_psi_del_phi_bary, ~, bin_psi] = histcounts(r_bary_dg_ci_final, psi_edges);

del_bary_phi_psi_mean = nan(size(hs_psi_del_phi_bary)); % mean \del\phi_B
del_bary_phi_psi_err = nan(size(hs_psi_del_phi_bary)); % std err of \del\phi_B
del_bary_phi_datapoints = nan(size(hs_psi_del_phi_bary));

psi_edges = psi_edges(1:end-1) + (psi_edges(2) - psi_edges(1))/2; %\Psi edges
psi_edges_fit = (psi_edges*180)/pi;

fig_psi_vs_del_phi = figure('Position', [300 300 1500 1500]);
xrange = [-185 185];

subplot(2,2,1)

del_phi_b_violin = nan(length(del_grp_phi_final),length(hs_psi_del_phi_bary));

for i = 1:length(hs_psi_del_phi_bary)

    id_temp = bin_psi == i;
    del_phi_b = del_grp_phi_final(id_temp'); % \del\phi_B for the corresponding \psi_BD

    mean_del_phi = mean(del_phi_b); % mean change in \phi
    err_del_phi = std(del_phi_b)/sqrt(length(del_phi_b)); % % std err change in \phi
    mean_del_phi = (mean_del_phi*180)/pi; % convert it to degree
    err_del_phi = (err_del_phi*180)/pi; % convert it to degree
    del_bary_phi_datapoints(i) = length(del_phi_b);
    del_bary_phi_psi_mean(i) = mean_del_phi;
    del_bary_phi_psi_err(i) = err_del_phi;

    del_phi_b_violin(1:length(del_phi_b),i) = (del_phi_b*180)/pi;

    scatter(ones(size(del_phi_b))*psi_edges_fit(i), (del_phi_b*180)/pi, 5, 'MarkerEdgeColor',[0.8 .8 .8],...
              'MarkerFaceColor',[0.8 .8 .8])
    hold on
    scatter(psi_edges_fit(i), mean_del_phi, 100, 'd', 'MarkerEdgeColor', 'r',...
              'MarkerFaceColor', 'r')
    hold on
    scatter(psi_edges_fit(i), (median(del_phi_b)*180)/pi, 100, 's', 'MarkerEdgeColor','m',...
              'MarkerFaceColor', 'm')

end

box('on')

set(gca, 'XLim', xrange, 'XTick', -180:30:180, 'YLim', [-185 185], 'YTick', -180:90:180, ...
    'FontSize', 14, 'LineWidth', 1)
xl = xlabel('Viewing angle of barycenter $\psi_{\rm BD} (^{o})$', 'Interpreter', 'latex', ...
    'FontName', 'Helvetica', 'FontSize', 20);
xl.Position(2) = xl.Position(2) - abs(xl.Position(2)*0.15);
yl = ylabel('Turing rate $\delta\phi_{\rm B}$', 'Interpreter', 'latex', ...
    'FontName', 'Helvetica', 'FontSize', 20);
tl = title('a', 'FontSize', 20, 'FontName', 'Helvetica', 'Color', 'k');
tl.Position(1) = yl.Position(1) - 10;

subplot(2,2,2)

req_dp = [1 2 3 4 10 11 12 13];
% req_dp_el = 5:9;
% del_phi_b_violin(:,req_dp_el);
% req_dp = 1:length(ci_edges_fit);

errorbar(psi_edges_fit(req_dp), del_bary_phi_psi_mean(req_dp), del_bary_phi_psi_err(req_dp), ...
    'o', 'color', [189 0 43]/sum([189 0 43]),...
    'linewidth', 1.5, 'MarkerSize', 10, 'MarkerEdgeColor', ...
    [189 0 43]/sum([189 0 43]),'MarkerFaceColor', [189 0 43]/sum([189 0 43]))

set(gca, 'XLim', xrange, 'XTick', -180:30:180, 'YLim', [-100 100], 'YTick', -100:40:100, ...
    'FontSize', 14, 'LineWidth', 1)
xl = xlabel('Viewing angle of barycenter $\psi_{\rm BD}$ ($^{o}$)', 'Interpreter', 'latex', ...
    'FontName', 'Helvetica', 'FontSize', 20);
% xl.Position(2) = xl.Position(2) - abs(xl.Position(2)*0.15);
yl = ylabel(['Mean turning rate $\delta\phi_{\rm B}$ $(^{o}/$', num2str(del_t*dt), '$\rm \,s )$'], 'Interpreter', 'latex', ...
    'FontName', 'Helvetica', 'FontSize', 20);
tl = title('b', 'FontSize', 20, 'FontName', 'Helvetica', 'Color', 'k');
tl.Position(1) = yl.Position(1) - 10;

%% Turning angle of dog

psi_edges = -pi-pi/12:pi/6:pi+pi/12; % \psi_DB edges

[hs_psi_del_phi_dog, ~, bin_psi] = histcounts(r_dg_bary_ci_final, psi_edges);

del_dog_phi_psi_mean = nan(size(hs_psi_del_phi_dog)); % mean \del\phi_D
del_dog_phi_psi_err = nan(size(hs_psi_del_phi_dog)); % std err \del\phi_D
del_dog_phi_datapoints = nan(size(hs_psi_del_phi_dog));

psi_edges = psi_edges(1:end-1) + (psi_edges(2) - psi_edges(1))/2; %\Psi edges
psi_edges_fit = (psi_edges*180)/pi;

subplot(2,2,3)

del_phi_d_violin = nan(length(del_dog_phi_final),length(hs_psi_del_phi_dog));

for i = 1:length(hs_psi_del_phi_dog)

    id_temp = bin_psi == i;
    del_phi_d = del_dog_phi_final(id_temp');

    mean_del_phi = mean(del_phi_d); % \del\phi_D for the corresponding \psi_DB
    err_del_phi = std(del_phi_d)/sqrt(length(del_phi_d)); % % std err change in \phi
    mean_del_phi = (mean_del_phi*180)/pi; % convert it to degree
    err_del_phi = (err_del_phi*180)/pi; % convert it to degree
    del_dog_phi_datapoints(i) = length(del_phi_d);
    del_dog_phi_psi_mean(i) = mean_del_phi;
    del_dog_phi_psi_err(i) = err_del_phi;

    del_phi_d_violin(1:length(del_phi_d),i) = (del_phi_d*180)/pi;

    scatter(ones(size(del_phi_d))*psi_edges_fit(i), (del_phi_d*180)/pi, 5, 'MarkerEdgeColor',[0.8 .8 .8],...
              'MarkerFaceColor',[0.8 .8 .8])
    hold on
    scatter(psi_edges_fit(i), mean_del_phi, 100, 'd', 'MarkerEdgeColor', 'r',...
              'MarkerFaceColor', 'r')
    hold on
    scatter(psi_edges_fit(i), (median(del_phi_d)*180)/pi, 100, 's', 'MarkerEdgeColor','m',...
              'MarkerFaceColor', 'm')

end

box('on')

set(gca, 'XLim', xrange, 'XTick', -180:30:180, 'YLim', [-185 185], 'YTick', -180:90:180, ...
    'FontSize', 14, 'LineWidth', 1)
xl = xlabel('Viewing angle of dog $\psi_{\rm DB} (^{o})$', 'Interpreter', 'latex', ...
    'FontName', 'Helvetica', 'FontSize', 20);
xl.Position(2) = xl.Position(2) - abs(xl.Position(2)*0.15);
yl = ylabel('Turing rate $\delta\phi_{\rm D}$', 'Interpreter', 'latex', ...
    'FontName', 'Helvetica', 'FontSize', 20);
tl = title('c', 'FontSize', 20, 'FontName', 'Helvetica', 'Color', 'k');
tl.Position(1) = yl.Position(1) - 10;

subplot(2,2,4)

req_dp = 3:11;
% req_dp = 1:length(ci_edges_fit);

errorbar(psi_edges_fit(req_dp), del_dog_phi_psi_mean(req_dp), del_dog_phi_psi_err(req_dp), 'o', 'color', [55 97 173]/sum([55 97 173]),...
    'linewidth', 1.5, 'MarkerSize',10,'MarkerEdgeColor', [55 97 173]/sum([55 97 173]), 'MarkerFaceColor', [55 97 173]/sum([55 97 173]))

set(gca, 'XLim', xrange, 'XTick', -180:30:180, 'YLim', [-105 110], 'YTick', -100:40:100, ...
    'FontSize', 14, 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k')
xl = xlabel('Viewing angle of Dog $\psi_{\rm DB}$ ($^{o}$)', 'Interpreter', 'latex', ...
    'FontName', 'Helvetica', 'FontSize', 20, 'Color', 'k');
xl.Position(2) = xl.Position(2) - abs(xl.Position(2)*0.15);
yl = ylabel(['Mean turning rate $\delta\phi_{\rm D}$ $(^{o}/$', num2str(del_t*dt), '$\rm s )$'], ...
    'Interpreter', 'latex', 'FontName', 'Helvetica', 'FontSize', 20, 'Color', 'k');
tl = title('d', 'FontSize', 20, 'FontName', 'Helvetica', 'Color', 'k');
tl.Position(1) = yl.Position(1) - 10;


%% plotting Fig.5 of main text

fig_5 = figure('Position', [300 300 700 1200]);

fs_label = 22;
fs_gca = 22;
fs_title = 30;

subplot(2,1,1)
violinplot(del_phi_b_violin, string(round(psi_edges_fit)), ...
    'ViolinColor', ones(size(del_phi_b_violin,2),3).*[0.6 0.6 .9], 'MedianMarkerSize', 70, ...
    'ViolinAlpha', 0.05, 'MarkerSize', 5)
set(gca, 'YLim', [-180 180], 'YTick', -180:60:180, ...
    'FontSize', fs_gca, 'FontName', font_name, 'LineWidth', 1, ...
    'Xcolor', 'k', 'YColor', 'k')
xlabel('Viewing angle of Barycenter (^o)', ...
    'FontName', font_name, 'FontSize', fs_label, 'Color', 'k')
ylabel(['Turning angle (^o/', num2str(del_t*dt), ' s)'], ...
    'FontName', font_name, 'FontSize', fs_label, 'Color', 'k')
% tl = title('a', 'FontSize', fs_title, 'FontName', font_name, 'Color', 'k');
% tl.Position(1) = -1.4;

% fig_del_phi_d_violin = figure('Position', [300 300 700 500]);
subplot(2,1,2)
violinplot(del_phi_d_violin, string(round(psi_edges_fit)), 'ViolinColor', ...
    ones(size(del_phi_d_violin,2),3).*[0.9 0.6 0.6], 'MedianMarkerSize', 50, 'ViolinAlpha', 0.1, 'MarkerSize', 5)
set(gca, 'YLim', [-180 180], 'YTick', -180:60:180, ...
    'FontSize', fs_gca, 'FontName', font_name, 'LineWidth', 1, ...
    'Xcolor', 'k', 'YColor', 'k')
xlabel('Viewing angle of Dog (^o)', ...
    'FontName', font_name, 'FontSize', fs_label, 'Color', 'k')
ylabel(['Turning angle (^o/', num2str(del_t*dt), ' s)'], ...
    'FontName', font_name, 'FontSize', fs_label, 'Color', 'k')
% tl = title('b', 'FontSize', fs_title, 'FontName', font_name, 'Color', 'k');
% tl.Position(1) = -1.4;

exportgraphics(fig_5, 'fig_5.pdf', 'ContentType', 'vector')
exportgraphics(fig_5, 'figure_5.jpeg')