%% storing the data from .dat file to .mat file in a way that is easily handled in Matlab

close all
clear
clc

%%

% loading the file which has all the data (processed by Ramon)
all_dat = readmatrix('sheepR.dat'); 

no_evnts = 11; % no of events (or different trials)
no_ind = 16; % no of individuals (sheep + dog + shepherd)
evnt_end = zeros(11,1); % when do each event end

% Where do each trial end (1st column of the data set)

for i = 1:no_evnts
    evnt_end_temp = find(all_dat(:,1) == i, 1, 'last'); % find the last indices of evnt_end_temp
    if isempty(evnt_end_temp)
        evnt_end(i,1) = NaN;
    else
        evnt_end(i,1) = evnt_end_temp;
    end
end

%%
% 
% Individuals time-steps within an event

evnt_end = evnt_end(~isnan(evnt_end)); % removing events that are not in data
no_evnts = length(evnt_end); % Final no.of events
ind_data = zeros(no_evnts,no_ind); % for a given event, what is the last time point for ind i before we move to j

for k = 1:no_evnts

    for i = 1:no_ind
        ind_data(k,i) = find(all_dat(1:evnt_end(k),2) == i, 1, 'last');
    end

end

%%
% 
% Saving the data as I need: pos (position) [x y]

sheep_dat = struct();

for ev = 1:no_evnts

    if ev == 1
        st_evt = 1;
    else
        st_evt = evnt_end(ev - 1) + 1;
    end

    tm = ind_data(ev,2) - ind_data(ev,1); % length of data points for a given ind for a given event

    pos = zeros(no_ind, 2, tm); % storing position

    all_dat_evt_1_x = all_dat(st_evt:evnt_end(ev),4); % x coordinate
    all_dat_evt_1_y = all_dat(st_evt:evnt_end(ev),5); % y coordinate 
    st_tm = 1;

    for i = 1:no_ind
        pos(i,1,:) = all_dat_evt_1_x(st_tm:(i*tm));
        pos(i,2,:) = all_dat_evt_1_y(st_tm:(i*tm));
        st_tm = (tm*i) + 1;
    end

    sheep_dat.(['pos_ev_', num2str(ev)]) = pos; % storing pos for a given event

end

save('sheep_dat.mat', '-struct', 'sheep_dat')