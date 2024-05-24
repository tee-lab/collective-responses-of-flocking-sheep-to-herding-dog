%% Run this code to reproduce Fig.4a-b of the main text. 

close all
clear
clc

%% %% Loading data 

load('ext_drives_data.mat') 
no_evnts = length(events);

no_shp = no_ind - 2; % no.of sheep, i.e., no.of individuals - (dog + shepherd)
font_name = 'Arial';
font_size = 20;

%% calculating pdfs of \psi_BD, \psi_DB, \phi_BD

r_bary_dg_ci_final = []; % \Psi_{BD}
r_dg_bary_ci_final = []; % \Psi_{DB}    
r_phi_bary_dg_final = []; % \phi_{BD}
del_t = 1;

for ev = 1:no_evnts

    pos = eval(strcat('pos_ev_', num2str(ev))); % load position data
    vel = eval(strcat('vel_ev_',num2str(ev))); % load velocity

    pos_shp = pos(1:14,:,2:del_t:end); % position of sheep 
    vel_shp = vel(1:14,:,2:del_t:end); % vel of sheep
    pos_dog = squeeze(pos(15,:,2:del_t:end)); % pos of dog
    dog_vel = squeeze(vel(15,:,2:del_t:end)); % vel of dog
    tm = size(pos_shp,3);
    
    r_bary_dg_ci = nan(1,tm-1); % where is dog from barycentre (\psi_BD)
    r_dg_bary_ci = nan(1,tm-1); % where is barycenter from dogs prespective (\psi_DB)
    r_phi_bary_dg = nan(1,tm-1); % orientation between bary and dog (\phi_BD)

    for t = 1:(tm-1)

        r_bary_dg_temp = pos_dog(:,t).' - mean(pos_shp(:,:,t),1); % displacement vec between bary and dog
        r_bary_dg_temp = r_bary_dg_temp/vecnorm(r_bary_dg_temp,2,2); % normalising 
        r_dg_bary_temp = -r_bary_dg_temp; % displacement vec between dog and barycenter

        vel_bary_temp = mean(vel_shp(:,:,t),1); % vel of bary
        vel_bary_temp = vel_bary_temp/(vecnorm(vel_bary_temp,2,2)+eps); % normalised vel of bary
        vel_dog_temp = dog_vel(:,t); % vel of dog
        vel_dog_temp = vel_dog_temp';
        vel_dog_temp = vel_dog_temp/(vecnorm(vel_dog_temp,2,2)+eps); % normalising

        theta_bary_dg = dot(vel_bary_temp, vel_dog_temp, 2); % v_B \cdot v_D
        theta_bary_dg = acos(theta_bary_dg); % find angle btw v_B and v_D
        theta_bary_dg_sign = cross([vel_bary_temp 0], [vel_dog_temp 0],2); % take cross product to check the sign (right hand thumb rule)
        theta_bary_dg = theta_bary_dg*sign(theta_bary_dg_sign(1,3)); % \phi_BD = \phi_D - \phi_B
        r_phi_bary_dg(1,t) = theta_bary_dg; % \phi_BD

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
    
    % concatenating all the data
    r_bary_dg_ci_final = [r_bary_dg_ci_final r_bary_dg_ci];
    r_dg_bary_ci_final = [r_dg_bary_ci_final r_dg_bary_ci];
    r_phi_bary_dg_final = [r_phi_bary_dg_final r_phi_bary_dg];

end

edges = linspace(-180, 180, 50); % bin edges for angles
r_bary_dg_ci_final = (180*r_bary_dg_ci_final)/pi; % \psi_BD
r_dg_bary_ci_final = (180*r_dg_bary_ci_final)/pi; % \psi_DB
r_phi_bary_dg_final = (180*r_phi_bary_dg_final)/pi; % \phi_BD

% calculating pdfs for all the above angles
[prob_den_bary_dg_ci, edges] = histcounts(r_bary_dg_ci_final, edges, 'Normalization', 'pdf');
[prob_den_dg_bary_ci, edges] = histcounts(r_dg_bary_ci_final, edges, 'Normalization', 'pdf');
[prob_den_bary_dg_phi, edges] = histcounts(r_phi_bary_dg_final, edges, 'Normalization', 'pdf');
edges = edges(1:end-1) + (edges(2)-edges(1))/2;

%% plotting Fig.4a-b of main text

fig_4a = figure('Position', [300 300 700 1200]);
subplot(2,1,1)

plot(edges, prob_den_bary_dg_ci, '-', 'LineWidth', 3, 'Color', '#0047AB')
hold on 
plot(edges, prob_den_dg_bary_ci, '-', 'LineWidth', 3, 'Color', '#C70039')

set(gca, 'XLim', [-180 180], 'XTick', -180:90:180, 'YLim', [0 0.015], 'YTick', 0:0.005:0.015, ...
    'YTickLabel', (0:0.005:0.015)*10^3, 'FontSize', font_size, 'FontName', font_name, ...
    'LineWidth', 1, 'Xcolor', 'k', 'YColor', 'k')
xl = xlabel('Viewing angle \psi', 'FontName', font_name, 'FontSize', font_size);
yl = ylabel('PDF (x 10^{-3})', 'FontName', font_name, 'FontSize', font_size);
% tl = title('a', 'FontSize', font_size+10, 'FontName', font_name, 'Color', 'k');
% tl.Position(1) = xl.Position(1) - 210;

legend({'$\psi_{\rm BD}$', '$\psi_{\rm DB}$'}, 'Interpreter', 'latex', ...
    'FontSize', font_size, 'Location', 'northwest')
legend('boxoff')


subplot(2,1,2)
plot(edges, prob_den_bary_dg_phi, '-', 'LineWidth', 3, 'Color', '#0047AB')
set(gca, 'XLim', [-180 180], 'XTick', -180:90:180, 'YLim', [0 0.012], 'YTick', 0:0.003:0.012, ...
    'YTickLabel', (0:0.003:0.012)*10^3, 'FontSize', font_size, 'FontName', font_name, ...
    'LineWidth', 1, 'Xcolor', 'k', 'YColor', 'k')
xl = xlabel('Heading difference \phi', 'FontName', font_name, 'FontSize', font_size);
yl = ylabel('PDF (x 10^{-3})', 'FontName', font_name, 'FontSize', font_size);
% tl = title('b', 'FontSize', font_size+10, 'FontName', font_name, 'Color', 'k');
% tl.Position(1) = xl.Position(1) - 210;

legend({'$\phi_{\rm BD}$'}, 'Interpreter', 'latex', ...
    'FontSize', font_size, 'Location', 'northwest')
legend('boxoff')
