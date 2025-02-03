%% Here I calculate heading direction, velocity and acceleration from position data

% If position data starts from t = 1, velocity and heading are valid only from t =
% 2 and accelerations are valid only from t = 3. 

% Run pre_ana.m file before running this file. 

close all
clear
clc

%%

load('sheep_dat.mat') 

no_ind = 16; % no of individuals
no_evnts = 9; % no of events
dt = 0.1; % time btw 2 measurements (0.1 s)

% calculating heading angles

sheep_all_dat = struct(); 

for ev = 1:no_evnts

    pos = eval(strcat('pos_ev_', num2str(ev))); % load position for given event
    pos_x = squeeze(pos(:,1,:)); % x coordinate for event i
    pos_y = squeeze(pos(:,2,:)); % y coordinate for event i
    tm = length(pos);

    vel = zeros(no_ind,2,tm); % storing velocity
    phi = zeros(no_ind, tm); % storing calculated heading angle

    for t = 2:tm

        vel_x = (pos_x(:,t) - pos_x(:,(t-1)))/dt; % velocity along x -> vx(t) = (x(t) - x(t-1))/dt
        vel_y = (pos_y(:,t) - pos_y(:,(t-1)))/dt; % similarly for y
        phi_temp = atan2(vel_y, vel_x); % calculating heading angle
        ind_phi = find(phi_temp < 0); % checking if anything is < 0 and setting everything btw [0,2pi]
        phi_temp(ind_phi) = phi_temp(ind_phi) + 2*pi;

        vel(:,1,t) = vel_x;
        vel(:,2,t) = vel_y;
        phi(:,t) = phi_temp;

    end
    
    vel_x = squeeze(vel(:,1,:));
    vel_y = squeeze(vel(:,2,:));

    sheep_all_dat.(['pos_ev_', num2str(ev)]) = pos; % storing pos and heading angle for a given event
    sheep_all_dat.(['vel_ev_', num2str(ev)]) = vel; % similarly for velocity
    sheep_all_dat.(['phi_ev_', num2str(ev)]) = phi; % for calculated heading angle
end

sheep_all_dat.no_evnts = no_evnts;
sheep_all_dat.no_ind = no_ind;
sheep_all_dat.dt = dt;

save('sheep_all_dat.mat', '-struct', 'sheep_all_dat')