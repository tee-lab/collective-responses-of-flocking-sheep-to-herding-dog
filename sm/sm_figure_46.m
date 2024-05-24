%% Run this code to reproduce Figures 3-5 of the SM

close all
clear
clc

%% Loading data 

load('ext_drives_data.mat') % can load both 'ext_drives_data.mat' and 'drives_data.mat'
load('sheep_all_dat.mat')
    
no_sp = ceil(sqrt(length(events))); % grids in subplot
tm_delay = 40;
c_min = 0.5;
font_name = 'Arial';
font_size = 20;

%% Calculate cross-correlation between group polarisation and median sheep speed

sm_fig_5a = figure('Position', [300 300 1000 1000]);

ccf_max = []; % cc coefficient 
ccf_lag = []; % time lag at ccc

for ev = 1:length(events)

    evt = events(ev);
    ev_st = eval(strcat('ev_st_', num2str(evt))); % start time 
    ev_et = eval(strcat('ev_et_', num2str(evt))); % end time
    drvs = length(ev_st); % no.of drives in the event

    vel = eval(strcat('vel_ev_',num2str(evt))); % load velocity
    phi_temp = eval(strcat('phi_ev_',num2str(evt))); % load heading angles (phi)

    ccf_max_temp = nan(1,drvs); % store max value of correlation
    ccf_lag_temp = nan(1,drvs); % store lag at max correlation

    for dr = 1:drvs

        ev_st_dr = ev_st(dr); % drive start time
        ev_et_dr = ev_et(dr); % drive end time

        vel_sheep = vel(1:14,:,ev_st_dr:ev_et_dr); % velocity of sheep
        spd_sheep = squeeze(vecnorm(vel_sheep,2,2)); % sheep speed
        
        % calculating median sheep speed
        median_shp_spd = median(spd_sheep, 1);
        median_shp_spd = median_shp_spd(2:end).';
        mean_med_shp_spd = mean(median_shp_spd);
        median_shp_spd = median_shp_spd - mean_med_shp_spd;

        % calculating group polarisation
        mx = mean(cos(phi_temp(1:14,ev_st_dr:ev_et_dr)),1); % mx
        mx = mx.';
        mx = mx(2:end); % not taking the 1st value because its always zero
        my = mean(sin(phi_temp(1:14,ev_st_dr:ev_et_dr)),1); % my
        my = my.';
        my = my(2:end); % not taking the 1st value because its always zero
        m = sqrt(mx.^2 + my.^2);
        mean_m = mean(m);
        m = m - mean_m;

        % calculating cross-corr between polarisation and median sheep
        % speed
        [ccf_m_mss_temp, tlag] = xcorr(m, median_shp_spd, tm_delay, 'normalized');
        
        ccf_m_mss_temp_abs = abs(ccf_m_mss_temp);
        tlag_id = find(ccf_m_mss_temp_abs == max(ccf_m_mss_temp_abs));

        len_data = length(median_shp_spd);
        ccf_signi = 2/sqrt(len_data - abs(tlag(tlag_id)));

        if max(abs(ccf_m_mss_temp)) >= c_min
            ccf_max_temp(1,dr) = ccf_m_mss_temp(tlag_id);
            ccf_lag_temp(1,dr) = tlag(tlag_id);
            subplot(no_sp,no_sp,ev)
            plot(tlag*dt, ccf_m_mss_temp, '-', 'LineWidth', 1)
            hold on
        end

    end

    ccf_max = [ccf_max ccf_max_temp];
    ccf_lag = [ccf_lag ccf_lag_temp];
    
    set(gca, 'FontName', font_name, 'FontSize', font_size, ...
        'XLim', [min(tlag*dt) max(tlag*dt)])
    xlabel('\tau (s)', 'FontName', font_name, 'FontSize', font_size)
    ylabel('C(\tau)', 'FontName', font_name, 'FontSize', font_size)

end

exportgraphics(sm_fig_5a, 'sm_fig_5a.pdf', 'ContentType', 'vector')

% how many drives cross-corr is greater than c_min
x_dat = sum(~isnan(ccf_lag));

% Plotting scatter plots of C_max and corresponding time lag.

sm_fig_5b = figure('Position', [300 300 600 500]);

if x_dat > 0

    % subplot(1,2,2)
    % h = gca;
    % x_dat_tlag = 1 + 0.02*randn(x_dat,1);
    % ccf_lag_plt = ccf_lag(~isnan(ccf_lag));
    % 
    % scatter(x_dat_tlag, ccf_lag_plt*dt, 40, 'MarkerEdgeColor', 'k',  ...
    % 'MarkerFaceColor', 'k')
    % hold on
    % boxplot(ccf_lag_plt*dt, 'Positions', 1, 'Widths', 0.2)
    % ylim([-0.25 0.25]);
    % h.YTick = -0.2:0.1:0.2;
    % h.XTickLabel = [];
    % h.FontName = font_name;
    % h.FontSize = font_size;
    % ylabel('Lag at C_{max} (s)')
    
    % subplot(1,2,1)
    h = gca;
    x_dat_cmax = 2 + 0.02*randn(x_dat,1);
    ccf_max_plt = ccf_max(~isnan(ccf_max));

    scatter(x_dat_cmax, ccf_max_plt, 40, 'MarkerEdgeColor', 'k',  ...
    'MarkerFaceColor', 'k')
    hold on
    boxplot(ccf_max_plt, 'Positions', 2, 'Widths', 0.2)
    ylim([0.5 0.9]);
    h.YTick = 0.5:0.1:0.9;
    h.XTickLabel = [];
    h.FontName = font_name;
    h.FontSize = font_size;
    ylabel('C_{max}')

end

exportgraphics(sm_fig_5b, 'sm_fig_5b.pdf', 'ContentType', 'vector')

%% Calculating cross-correlation between dog speed and median speed

ccf_max = []; % cc coefficient 
ccf_lag = []; % time lag at ccc

sm_fig_6a = figure('Position', [300 300 1000 1000]);
plt = 1;
plt_dr = 0;

for ev = 1:length(events)

    evt = events(ev);
    ev_st = eval(strcat('ev_st_', num2str(evt))); % start time 
    ev_et = eval(strcat('ev_et_', num2str(evt))); % end time
    drvs = length(ev_st); % no.of drives in the event

    vel = eval(strcat('vel_ev_',num2str(evt))); % load velocity

    ccf_max_temp = nan(1,drvs); % storing C_max
    ccf_lag_temp = nan(1,drvs); % storting t_lag

    for dr = 1:drvs

        ev_st_dr = ev_st(dr); % drive start time
        ev_et_dr = ev_et(dr); % drive end time

        vel_sheep = vel(1:14,:,ev_st_dr:ev_et_dr); % velocity of sheep
        vel_dog = squeeze(vel(15,:,ev_st_dr:ev_et_dr)); % velocity of dog

        spd_sheep = squeeze(sqrt(vel_sheep(:,1,:).^2 + vel_sheep(:,2,:).^2)); % sheep speed
        spd_dog = sqrt(vel_dog(1,:).^2 + vel_dog(2,:).^2); % dog speed
        spd_dog = spd_dog(2:end);

        vel_grp = squeeze(mean(vel_sheep,1));
        spd_grp = vecnorm(vel_grp, 2, 1);
        spd_grp = spd_grp - mean(spd_grp);
        spd_grp = spd_grp(2:end);
        
        % calculating median dog speed
        median_shp_spd = median(spd_sheep, 1);
        median_shp_spd = median_shp_spd(2:end).';
        mean_med_shp_spd = mean(median_shp_spd);
        median_shp_spd = median_shp_spd - mean_med_shp_spd;
        
        % Calculating dog speed
        spd_dog = spd_dog - mean(spd_dog);

        % cross-corr between dog speed and median sheep speed
        % [ccf_dgs_mss_temp, tlag] = xcorr(spd_dog, median_shp_spd, tm_delay, 'normalized');
        [ccf_dgs_mss_temp, tlag] = xcorr(spd_dog, spd_grp, tm_delay, 'normalized');
        
        ccf_dgs_mss_temp_abs = abs(ccf_dgs_mss_temp);
        tlag_id = find(ccf_dgs_mss_temp_abs == max(ccf_dgs_mss_temp_abs));

        len_data = length(median_shp_spd);
        ccf_signi = 2/sqrt(len_data - abs(tlag(tlag_id)));
        
        if max(abs(ccf_dgs_mss_temp)) >= c_min  
            ccf_max_temp(1,dr) = ccf_dgs_mss_temp(tlag_id);
            ccf_lag_temp(1,dr) = tlag(tlag_id);
            subplot(3,2,plt)
            plot(tlag*dt, ccf_dgs_mss_temp, '-', 'LineWidth', 1)
            plt_dr = 1;
            hold on
        end

    end

    if plt_dr == 1
        plt = plt+1;
        plt_dr = 0;
    end

    ccf_max = [ccf_max ccf_max_temp];
    ccf_lag = [ccf_lag ccf_lag_temp];

    set(gca, 'FontName', font_name, 'FontSize', font_size, ...
        'XLim', [min(tlag*dt) max(tlag*dt)])
    xlabel('\tau (s)', 'FontName', font_name, 'FontSize', font_size)
    ylabel('C(\tau)', 'FontName', font_name, 'FontSize', font_size)

end

exportgraphics(sm_fig_6a, 'sm_fig_6a.pdf', 'ContentType', 'vector')

x_dat = sum(~isnan(ccf_lag));

sm_fig_6b = figure('Position', [300 300 600 500]);

if x_dat > 0

    % subplot(1,2,2)
    % h = gca;
    % 
    % x_dat_tlag = 1 + 0.02*randn(x_dat,1);
    % ccf_lag_plt = ccf_lag(~isnan(ccf_lag));
    % 
    % scatter(x_dat_tlag, ccf_lag_plt*dt, 40, 'MarkerEdgeColor', 'k',  ...
    % 'MarkerFaceColor', 'k')
    % hold on
    % boxplot(ccf_lag_plt*dt, 'Positions', 1, 'Widths', 0.2)
    % % ylim([-4.5 4.5])
    % ylabel('Lag at C_{max} (s)')
    % h.YAxis.TickValues = -1:1:2;
    % h.YAxis.TickLength = [0.007 0.007];
    % h.XAxis.TickValues = [];
    % h.FontName = font_name;
    % h.FontSize = font_size;
    
    % subplot(1,2,1)
    h = gca;

    x_dat_cmax = 2 + 0.02*randn(x_dat,1);
    ccf_max_plt = ccf_max(~isnan(ccf_max));

    scatter(x_dat_cmax, ccf_max_plt, 40, 'MarkerEdgeColor', 'k',  ...
    'MarkerFaceColor', 'k')
    hold on
    boxplot(ccf_max_plt, 'Positions', 2, 'Widths', 0.2)
    ylim([0.5 0.9])
    ylabel('C_{max}')
    h.YAxis.TickValues = 0.5:.1:0.9;
    h.YAxis.TickLength = [0.007 0.007];
    h.XAxis.TickValues = [];
    h.FontName = font_name;
    h.FontSize = font_size;

end

exportgraphics(sm_fig_6b, 'sm_fig_6b.pdf', 'ContentType', 'vector')
%% Calculate cross-correlation between dog velocity and grp velocity

ccf_max = []; % cc coefficient 
ccf_lag = []; % time lag at ccc

sm_fig_7a = figure('Position', [300 300 1000 1000]);

for ev = 1:length(events)

    evt = events(ev);
    ev_st = eval(strcat('ev_st_', num2str(evt))); % start time 
    ev_et = eval(strcat('ev_et_', num2str(evt))); % end time
    drvs = length(ev_st); % no.of drives in the event

    vel = eval(strcat('vel_ev_',num2str(evt))); % load velocity

    ccf_max_temp = nan(1,drvs); % storing c_max
    ccf_lag_temp = nan(1,drvs); % storing t_lag

    for dr = 1:drvs

        ev_st_dr = ev_st(dr); % drive start time
        ev_et_dr = ev_et(dr); % drive end time

        vel_sheep = vel(1:14,:,ev_st_dr:ev_et_dr); % velocity of sheep
        grp_vel = squeeze(mean(vel_sheep,1));
        grp_vel = grp_vel./sqrt(grp_vel(1,:).^2 + grp_vel(2,:).^2);
        vel_dog = squeeze(vel(15,:,ev_st_dr:ev_et_dr)); % velocity of dog
        vel_dog = vel_dog./sqrt(vel_dog(1,:).^2 + vel_dog(2,:).^2);
        
        % calculating cross-correlation between dog velocity and v_B
        [ccf_vdx_vgx_temp, ~] = xcorr(vel_dog(1,2:end), grp_vel(1,2:end), tm_delay, 'unbiased'); % x component
        [ccf_vdy_vgy_temp, tlag] = xcorr(vel_dog(2,2:end), grp_vel(2,2:end), tm_delay, 'unbiased'); % y component

        ccf_vel_temp = (ccf_vdx_vgx_temp + ccf_vdy_vgy_temp); % x + y
        
        ccf_vel_temp_abs = abs(ccf_vel_temp);
        tlag_id = find(ccf_vel_temp_abs == max(ccf_vel_temp_abs));

        ccf_signi = 0.5;
        
        if max(abs(ccf_vel_temp)) >= ccf_signi
            ccf_max_temp(1,dr) = ccf_vel_temp(tlag_id);
            ccf_lag_temp(1,dr) = tlag(tlag_id);
            subplot(no_sp,no_sp,ev)
            plot(tlag*dt, ccf_vel_temp, '-', 'LineWidth', 1)
            hold on
        end

    end

    ccf_max = [ccf_max ccf_max_temp];
    ccf_lag = [ccf_lag ccf_lag_temp];
    
     set(gca, 'FontName', font_name, 'FontSize', font_size, ...
        'XLim', [min(tlag*dt) max(tlag*dt)])
    xlabel('\tau (s)', 'FontName', font_name, 'FontSize', font_size)
    ylabel('C(\tau)', 'FontName', font_name, 'FontSize', font_size)

end

exportgraphics(sm_fig_7a, 'sm_fig_7a.pdf', 'ContentType', 'vector')

sm_fig_7b = figure('Position', [300 300 600 500]);

x_dat = sum(~isnan(ccf_lag));

if x_dat > 0

    % subplot(1,2,2)
    % h = gca;
    % 
    % x_dat_tlag = 1 + 0.02*randn(x_dat,1);
    % ccf_lag_plt = ccf_lag(~isnan(ccf_lag));
    % 
    % scatter(x_dat_tlag, ccf_lag_plt*dt, 40, 'MarkerEdgeColor', 'k',  ...
    % 'MarkerFaceColor', 'k')
    % hold on
    % boxplot(ccf_lag_plt*dt, 'Positions', 1, 'Widths', 0.2)
    % ylim([-4.9 4.2])
    % ylabel('Lag at C_{max} (s)')
    % h.YAxis.TickValues = -4:2:4;
    % h.YAxis.TickLength = [0.007 0.007];
    % h.XAxis.TickValues = [];
    % h.FontName = font_name;
    % h.FontSize = font_size;
    
    % subplot(1,2,1)
    h = gca;

    x_dat_cmax = 2 + 0.02*randn(x_dat,1);
    ccf_max_plt = ccf_max(~isnan(ccf_max));

    scatter(x_dat_cmax, ccf_max_plt, 40, 'MarkerEdgeColor', 'k',  ...
    'MarkerFaceColor', 'k')
    hold on
    boxplot(ccf_max_plt, 'Positions', 2, 'Widths', 0.2)
    ylim([0.7 1])
    ylabel('C_{max}')
    h.YAxis.TickValues = 0.7:.1:1;
    h.YAxis.TickLength = [0.007 0.007];
    h.XAxis.TickValues = [];
    h.FontName = font_name;
    h.FontSize = font_size;

end

exportgraphics(sm_fig_7b, 'sm_fig_7b.pdf', 'ContentType', 'vector')