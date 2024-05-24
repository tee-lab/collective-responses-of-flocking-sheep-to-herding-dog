%% Collecting data for visually identified drives

close all
clear
clc

%% Loading data

load('sheep_all_dat.mat')
drive_data = struct();

% drives for a given event
ev_st_4 = [50 130 214]/dt;
ev_et_4 = [77 170 240]/dt;

ev_st_5 = [0.1 48 110]/dt;
ev_et_5 = [24 68 122]/dt;

ev_st_8 = [0.1 55 170 332 438]/dt;
ev_et_8 = [33 90 220 387 530]/dt;

ev_st_9 = [0.1 99]/dt;
ev_et_9 = [32 154]/dt;

% events in which drives were observed
events = [4 5 8 9];

for evt = 1:length(events)

    ev = events(evt); % event of the drive

    ev_st_temp = eval(strcat('ev_st_', num2str(ev))); % list of start of drive
    ev_et_temp = eval(strcat('ev_et_', num2str(ev))); % list of end of drive

    drive_data.(['ev_st_', num2str(ev)]) = ev_st_temp; % storing start of drive 
    drive_data.(['ev_et_', num2str(ev)]) = ev_et_temp; % storing end times of drive

    phi = eval(strcat('phi_ev_',num2str(ev))); % load heading angles (phi)
    vel = eval(strcat('vel_ev_',num2str(ev))); % load velocity
    pos = eval(strcat('pos_ev_', num2str(ev))); % load position

    phi_drive = [];
    vel_drive_x = [];
    vel_drive_y = [];
    pos_drive_x = [];
    pos_drive_y = [];
    a_drive_x = [];
    a_drive_y = [];
    
    % collecting data over all drives (and only drives)
    for tdr = 1:length(ev_st_temp)

        st = ev_st_temp(tdr);
        et = ev_et_temp(tdr);

        phi_drive = [phi_drive phi(:,st:et)];
        vel_drive_x = [vel_drive_x squeeze(vel(:,1,st:et))];
        vel_drive_y = [vel_drive_y squeeze(vel(:,2,st:et))];
        pos_drive_x = [pos_drive_x squeeze(pos(:,1,st:et))];
        pos_drive_y = [pos_drive_y squeeze(pos(:,2,st:et))];
        a_drive_x = [a_drive_x squeeze(vel(:,1,st:et))];
        a_drive_y = [a_drive_y squeeze(vel(:,2,st:et))];

    end

    tt = length(phi_drive);
    vel_drive = zeros(no_ind, 2, tt);
    vel_drive(:,1,:) = vel_drive_x;
    vel_drive(:,2,:) = vel_drive_y;

    pos_drive = zeros(no_ind, 2, tt);
    pos_drive(:,1,:) = pos_drive_x;
    pos_drive(:,2,:) = pos_drive_y;

    acl_drive = zeros(no_ind, 2, tt);
    acl_drive(:,1,:) = a_drive_x;
    acl_drive(:,2,:) = a_drive_y;

    drive_data.(['phi_ev_', num2str(evt)]) = phi_drive;
    drive_data.(['pos_ev_', num2str(evt)]) = pos_drive;
    drive_data.(['vel_ev_', num2str(evt)]) = vel_drive;
    drive_data.(['acl_ev_', num2str(evt)]) = acl_drive;

end

drive_data.events = events; 
drive_data.no_ind = no_ind;
drive_data.dt = dt;

save('drives_data.mat', '-struct', 'drive_data')