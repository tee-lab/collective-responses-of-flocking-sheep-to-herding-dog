%% Collecting data for visually identified drives

close all
clear
clc

%% Loading data

load('sheep_all_dat.mat') % load data
drive_data = struct();

% identified drive times from each event (or, trail)

ev_st_1 = [140 453]/dt; % Start time of drives
ev_et_1 = [160 465]/dt; % End time of drives
len_drvs_evt_1 = ev_et_1 - ev_st_1; % length of drives

ev_st_2 = [151 196 258 700]/dt;
ev_et_2 = [177 232 290 730]/dt;
len_drvs_evt_2 = ev_et_2 - ev_st_2;

ev_st_3 = 163/dt;
ev_et_3 = 193/dt;
len_drvs_evt_3 = ev_et_3 - ev_st_3;

ev_st_4 = [23 50 130 191 210 272 293]/dt;
ev_et_4 = [35 77 170 198 240 286 298]/dt;
len_drvs_evt_4 = ev_et_4 - ev_st_4;

ev_st_5 = [0.1 48 110 207 283]/dt;
ev_et_5 = [24 68 128 227 288]/dt;
len_drvs_evt_5 = ev_et_5 - ev_st_5;

ev_st_6 = 81/dt;
ev_et_6 = 89/dt;
len_drvs_evt_6 = ev_et_6 - ev_st_6;

ev_st_7 = [];
ev_et_7 = [];

ev_st_8 = [0.1 55 170 332 438 585]/dt;
ev_et_8 = [33 90 220 387 530 593]/dt;
len_drvs_evt_8 = ev_et_8 - ev_st_8;

ev_st_9 = [0.1 99 303]/dt;
ev_et_9 = [32 210 324]/dt;
len_drvs_evt_9 = ev_et_9 - ev_st_9;

len_drvs = [len_drvs_evt_1 len_drvs_evt_2 len_drvs_evt_3 len_drvs_evt_4 len_drvs_evt_5 len_drvs_evt_6 len_drvs_evt_8 len_drvs_evt_9];

events = [1 2 3 4 5 6 8 9];

for evt = 1:length(events)

    ev = events(evt); % Event or trial number

    ev_st_temp = eval(strcat('ev_st_', num2str(ev))); % start times of drives
    ev_et_temp = eval(strcat('ev_et_', num2str(ev))); % end times of drives

    drive_data.(['ev_st_', num2str(ev)]) = ev_st_temp; % store start times
    drive_data.(['ev_et_', num2str(ev)]) = ev_et_temp; % store end times 

    phi = eval(strcat('phi_ev_',num2str(ev))); % load heading angles (phi)
    vel = eval(strcat('vel_ev_',num2str(ev))); % load velocity
    pos = eval(strcat('pos_ev_', num2str(ev))); % load position

    % storing heading angles, vel (x and y), pos (x and y) for all drives
    % in a given event. 
    phi_drive = []; 
    vel_drive_x = [];
    vel_drive_y = [];
    pos_drive_x = [];
    pos_drive_y = [];

    for tdr = 1:length(ev_st_temp)

        st = ev_st_temp(tdr);
        et = ev_et_temp(tdr);

        phi_drive = [phi_drive phi(:,st:et)];
        vel_drive_x = [vel_drive_x squeeze(vel(:,1,st:et))];
        vel_drive_y = [vel_drive_y squeeze(vel(:,2,st:et))];
        pos_drive_x = [pos_drive_x squeeze(pos(:,1,st:et))];
        pos_drive_y = [pos_drive_y squeeze(pos(:,2,st:et))];

    end

    tt = length(phi_drive);
    vel_drive = zeros(no_ind, 2, tt);
    vel_drive(:,1,:) = vel_drive_x;
    vel_drive(:,2,:) = vel_drive_y;

    pos_drive = zeros(no_ind, 2, tt);
    pos_drive(:,1,:) = pos_drive_x;
    pos_drive(:,2,:) = pos_drive_y;

    drive_data.(['phi_ev_', num2str(evt)]) = phi_drive;
    drive_data.(['pos_ev_', num2str(evt)]) = pos_drive;
    drive_data.(['vel_ev_', num2str(evt)]) = vel_drive;

end

drive_data.events = events; 
drive_data.no_ind = no_ind;
drive_data.dt = dt;
drive_data.len_drvs = len_drvs;

save('ext_drives_data.mat', '-struct', 'drive_data')