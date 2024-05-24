function [pos_s_dat, pos_d_dat, vel_s_dat, vel_d_dat, spd_d_dat, collect_t, ...
    drive_t, force_slow_t] = herding_model(no_shp, box_length, ...
    rad_rep_s, rad_rep_dog, K_atr, k_atr, k_alg, vs, v_dog, h, rho_a, rho_d, e, c, alg_str, ...
    f_n, pd, pc, n_iter)

% initialisation

theta_pos = 2*pi*rand();
str_side = box_length*[cos(theta_pos) sin(theta_pos)];
pos_s = str_side - 3*rad_rep_s*rand(no_shp,2); % put sheep on top right corner of box
% pos_d = [rand() box_length - rand()]; % place dog on top left corner of box 
pos_d = str_side - 3*rad_rep_dog*rand(1,2);

theta_s = 2*pi*rand(no_shp,1);
theta_d = 2*pi*rand();
vel_s = [cos(theta_s) sin(theta_s)]; % initial sheep velocity
vel_d = [cos(theta_d) sin(theta_d)]; % ini dog velocity
spd_s = ones(no_shp,1)*vs;

% saving data 

pos_s_dat = nan(no_shp,2,n_iter); % position of sheep
vel_s_dat = nan(no_shp,2,n_iter); % velocity of sheep
pos_d_dat = nan(n_iter,2); % position of dog
vel_d_dat = nan(n_iter,2); % velocity of dog
spd_d_dat = nan(n_iter,1); % speed of dog
collect_t = nan(n_iter,1); % no of times dog collect
drive_t = nan(n_iter,1); % no of times dog drives
force_slow_t = nan(n_iter,1); % no of times dog is within r_a

% storing initial position and orientation
pos_s_dat(:,:,1) = pos_s;
vel_s_dat(:,:,1) = vel_s;
pos_d_dat(1,:) = pos_d;
vel_d_dat(1,:) = vel_d;
spd_d_dat(1) = v_dog;

%%

% position and orientation at time t-1.

pos_s_t_1 = pos_s;
pos_d_t_1 = pos_d;
vel_s_t_1 = vel_s;
vel_d_t_1 = vel_d;

for t = 2:n_iter

    % movement of sheep
    for i = 1:no_shp

        r_shp_dg = pos_d_t_1 - pos_s_t_1(i,:);
        dist_rsd = vecnorm(r_shp_dg);
        r_shp_dg = r_shp_dg/dist_rsd;

        if dist_rsd > rad_rep_dog % if the sheep beyond interaction radius with dog

            r_ij = pos_s_t_1 - repmat(pos_s_t_1(i,:), no_shp, 1);
            mag_rij = vecnorm(r_ij, 2, 2); % distance between sheep in the herd. 

            rep_j = find(mag_rij < rad_rep_s);

            if length(rep_j) > 1 % if there are any individuals within the repulsion zone
                shp_i = find(rep_j == i);
                rep_j(shp_i) = [];
                r_ij_rep = r_ij(rep_j,:)./mag_rij(rep_j);
                r_ij_rep = sum(r_ij_rep, 1);
                r_ij_rep = r_ij_rep/vecnorm(r_ij_rep);
                r_ij_rep = -r_ij_rep;

                vel_next = h*vel_s_t_1(i,:) + rho_a*r_ij_rep;
                vel_next = vel_next/vecnorm(vel_next);

                pos_s(i,:) = pos_s(i,:) + vs*vel_next;
                vel_s_dat(i,:,t) = vel_next;
            else
                pos_s(i,:) = pos_s(i,:);
                vel_s_dat(i,:,t) = vel_s_t_1(i,:);
            end

        else % when dog is seen by the sheep

            % repulsion between sheep
            r_ij = pos_s_t_1 - repmat(pos_s_t_1(i,:), no_shp, 1);
            mag_rij = vecnorm(r_ij, 2, 2);

            rep_j = find(mag_rij < rad_rep_s);
            is_err = 0;
            if length(rep_j) > 1
                shp_i = find(rep_j == i);
                rep_j(shp_i) = [];
                r_ij_rep = r_ij(rep_j,:)./mag_rij(rep_j);
                r_ij_rep = sum(r_ij_rep, 1);
                r_ij_rep = r_ij_rep/vecnorm(r_ij_rep);
                r_ij_rep = -r_ij_rep;
                is_err = 1;
            end

            % repulsion from dog
            r_shp_dg = -r_shp_dg;

            % attraction towards LCM
            [~, lcm_j] = sort(mag_rij, 'ascend');
            lcm_j = lcm_j(2:(K_atr+1));
            lcm_j = lcm_j(randperm(K_atr,k_atr));
            r_atr = r_ij(lcm_j,:)./(mag_rij(lcm_j)+eps);
            r_atr = sum(r_atr, 1);
            r_atr = r_atr/vecnorm(r_atr);

            % alignment
            l_alg = lcm_j(randperm(length(lcm_j),k_alg));
            r_alg = vel_s_t_1(l_alg,:);
            r_alg = sum(r_alg);
            r_alg = r_alg/vecnorm(r_alg);
%             s_alg = mean(spd_s_t_1(l_alg)); % speed is average of near neighbours

            % error in copying
            theta_error = 2*pi*rand();
            r_err = [cos(theta_error) sin(theta_error)];

            % resultant velocity
            if is_err == 1
                vel_next = h*vel_s_t_1(i,:) + rho_a*r_ij_rep +...
                    rho_d*r_shp_dg + c*r_atr + e*r_err + alg_str*r_alg;
            else
                vel_next = h*vel_s_t_1(i,:) + rho_d*r_shp_dg + c*r_atr + e*r_err + alg_str*r_alg;
            end

            vel_next = vel_next/vecnorm(vel_next);
            pos_s(i,:) = pos_s(i,:) + vs*vel_next;
            vel_s_dat(i,:,t) = vel_next;
            spd_s(i) = vs;

        end

    end

    % Movement of dog

    r_dg_shp = pos_s_t_1 - repmat(pos_d_t_1, no_shp, 1);
    dist_rds = vecnorm(r_dg_shp,2,2);

    % if the dog is within r_a distance from sheep
    if min(dist_rds) <= rad_rep_s

        pos_d = pos_d + 0.05*vel_d_t_1;
        vel_d_dat(t,:) = vel_d_t_1;
        spd_d_dat(t) = 0.05;
        force_slow_t(t) = 1;

    else

        grp_centre = mean(pos_s_t_1, 1);
        r_gcm_i = pos_s_t_1 - repmat(grp_centre, no_shp, 1);
        dist_gcm_i = vecnorm(r_gcm_i, 2, 2);
        
        % if group is not cohesive
        if max(dist_gcm_i) >= f_n

            s_p = find(dist_gcm_i == max(dist_gcm_i));
            d_behind = dist_gcm_i(s_p) + pc;
            rc = d_behind*(r_gcm_i(s_p,:)/dist_gcm_i(s_p));
            rc = grp_centre + rc;
            rdc = rc - pos_d_t_1;
            rdc = rdc/vecnorm(rdc);

            % error in copying
            theta_error = 2*pi*rand();
            r_err = [cos(theta_error) sin(theta_error)];

            vel_next = rdc + e*r_err;
            vel_next = vel_next/vecnorm(vel_next);

            pos_d = pos_d + v_dog*vel_next;
            collect_t(t) = 1;
        
        % if group is cohesive
        else 
            
            d_behind = vecnorm(grp_centre) + pd;
            r_drive = d_behind*(grp_centre/vecnorm(grp_centre));
            r_drive_orient = r_drive - pos_d_t_1;
            r_drive_orient = r_drive_orient/vecnorm(r_drive_orient);

            % error in copying
            theta_error = 2*pi*rand();
            r_err = [cos(theta_error) sin(theta_error)];

            vel_next = r_drive_orient + e*r_err;
            vel_next = vel_next/vecnorm(vel_next);

            pos_d = pos_d + v_dog*vel_next;
            drive_t(t) = 1;

        end

        vel_d_t_1 = vel_next;
        vel_d_dat(t,:) = vel_next;
        spd_d_dat(t) = v_dog;

    end

    pos_s_dat(:,:,t) = pos_s;
    pos_s_t_1 = pos_s;
    vel_s_t_1 = vel_s_dat(:,:,t);

    pos_d_dat(t,:) = pos_d;
    pos_d_t_1 = pos_d;

end