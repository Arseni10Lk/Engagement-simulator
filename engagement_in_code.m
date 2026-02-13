clc
clear all
close all

% dbstop if naninf


% Visualisation data
engagement_plots = 'on';
full_scale_engagement_viz = 'on';
zoomed_engagement_viz = 'on';

vid_file = 'on';

% initial conditions
PN_type = 'Pure';
at = [0, 0, 0]; % x y z

HE_pitch = -20*pi/180;
HE_yaw = -20*pi/180;

N = 3;
time_final = 12;
t_step = 0.01;

target_pitch = 20*pi/180;
target_yaw = 0;

Rt = [40000, 10000, 0]; % x y z [m]
Rp = [0    , 10000, 0]; % x y z [m]

VP_mag = 4000; % m/s
Vt_mag = 1000; % m/s

% calculation of dependent inital conditions
Vt = zeros(1, 3);
Vt(1) = -Vt_mag*cos(target_pitch)*cos(target_yaw);
Vt(2) = Vt_mag*sin(target_pitch);
Vt(3) = Vt_mag*cos(target_pitch)*sin(target_yaw);

RTM = Rt - Rp;

lambda_pitch = atan2(RTM(2), RTM(1)); % rad
lambda_yaw = atan2(RTM(3), RTM(1)); % rad

Lead_angle_pitch = asin(Vt_mag*sin(target_pitch + lambda_pitch)/VP_mag); % rad
Lead_angle_yaw = asin(Vt_mag*cos(target_pitch)*sin(target_yaw + lambda_yaw)/VP_mag/cos(Lead_angle_pitch)); % rad

Vpx = VP_mag*cos(Lead_angle_pitch + lambda_pitch + HE_pitch)*cos(Lead_angle_yaw + lambda_yaw + HE_yaw); % m/s 
Vpy = VP_mag*sin(Lead_angle_pitch + lambda_pitch + HE_pitch); % m/s
Vpz = VP_mag*cos(Lead_angle_pitch + lambda_pitch + HE_pitch)*sin(Lead_angle_yaw + lambda_yaw + HE_yaw); % m/s
Vp = [Vpx, Vpy, Vpz];

% Vector with Initial conditions

y0 = [Rt, Rp, Vt, Vp]';

% application of ProNav laws
time = 0:t_step:time_final-t_step; 
length_time = length(time);

Con_vs_t = zeros(length(y0), length_time); 

% Integration with RK4
Actual_state = y0;
Con_vs_t(:, 1) = Actual_state;

prev_range = norm(y0(1:3) - y0(4:6));

for j = 1 : length_time-1

    s1 = nlinpronav_sim(Actual_state, N, at, PN_type);
    s2 = nlinpronav_sim(Actual_state + t_step*s1/2, N, at, PN_type);
    s3 = nlinpronav_sim(Actual_state + t_step*s2/2, N, at, PN_type);
    s4 = nlinpronav_sim(Actual_state + t_step*s3, N, at, PN_type);

    Actual_state = Actual_state + t_step*(s1+2*s2+2*s3+s4)/6;
    
    Rt_curr = Actual_state(1:3);
    Rp_curr = Actual_state(4:6);
    curr_range = norm(Rt_curr - Rp_curr);

    if curr_range > prev_range
        % Trim the empty zeros from the pre-allocated array
        Con_vs_t = Con_vs_t(:, 1:j);
        disp(['Intercept occurred at t = ' num2str(j*t_step) ' sec']);
        break;
    end

    Con_vs_t(:, j + 1) = Actual_state;

    prev_range = curr_range;
end


Actual_state_history = Con_vs_t';

% Processing simulation results

sel_RT = 1:3;
sel_RM = 4:6;
sel_VT = 7:9;
sel_VM = 10:12;

VP_mag = vecnorm(Actual_state_history(:, sel_VM), 2, 2);
Vt_mag = vecnorm(Actual_state_history(:, sel_VT), 2, 2);

RTM = Actual_state_history(:, sel_RT)-Actual_state_history(:, sel_RM);
VTM = Actual_state_history(:, sel_VT)-Actual_state_history(:, sel_VM);

LOS_rate = cross(RTM, VTM, 2)./vecnorm(RTM, 2, 2).^2;

VC_hist = -dot(RTM, VTM, 2)./vecnorm(RTM, 2, 2);

if strcmp(PN_type, 'True')
    ref_vecs = RTM;
    vel_term = VC_hist;
elseif strcmp(PN_type, 'Pure')
    ref_vecs = Actual_state_history(:, sel_VM);
    vel_term = vecnorm(ref_vecs, 2, 2); % Missile Speed
else
    error("Wrong type");
end

% Normalize Reference Vectors (Nx3)
u_fwd_hist = ref_vecs ./ vecnorm(ref_vecs, 2, 2);
global_up = [0, 1, 0]; % 1x3 Row for expansion

% Calculate Right Vector (Nx3)
u_right_hist = cross(u_fwd_hist, repmat(global_up, length(u_fwd_hist), 1), 2);
u_right_hist = u_right_hist ./ vecnorm(u_right_hist, 2, 2);

% Calculate Local Up Vector (Nx3)
u_local_up_hist = cross(u_right_hist, u_fwd_hist, 2);
u_local_up_hist = u_local_up_hist ./ vecnorm(u_local_up_hist, 2, 2);

% 5. Project LOS Rate onto the Local Axes
%    Now we can correctly dot the Nx3 rate with the Nx3 basis vectors
lambda_dot_pitch = dot(LOS_rate, u_right_hist, 2);
lambda_dot_yaw   = dot(LOS_rate, u_local_up_hist, 2);

aM = N * abs(vel_term) .* vecnorm(LOS_rate, 2, 2);

sel_Rtx = 1; sel_Rty = 2; sel_Rtz = 3;
sel_Rpx = 4; sel_Rpy = 5; sel_Rpz = 6;
sel_Vtx = 7; sel_Vty = 8; sel_Vtz = 9;
sel_Vpx = 10; sel_Vpy = 11; sel_Vpz = 12;

time = time(1:length(aM));


dt_index = .1/t_step;
% Find minimum distance index (3D)
[min_dist, miss_index] = min(vecnorm(RTM, 2, 2));

if strcmp(engagement_plots, 'on')

    figure(1)
    
    % 1. Plot Missile Trajectory (Blue)
    % Note: Using (:, sel_Rpz) for Visual Y, and (:, sel_Rpy) for Visual Z (Altitude)
    plot3(Actual_state_history(1:miss_index, sel_Rpx), ...
          Actual_state_history(1:miss_index, sel_Rpz), ...
          Actual_state_history(1:miss_index, sel_Rpy), ...
          'b', 'LineWidth', 2)
    hold on 
    
    % 2. Plot Target Trajectory (Red Dashed)
    plot3(Actual_state_history(1:miss_index, sel_Rtx), ...
          Actual_state_history(1:miss_index, sel_Rtz), ...
          Actual_state_history(1:miss_index, sel_Rty), ...
          'r--', 'LineWidth', 2)
          
    % 3. Plot Start Points (Circles)
    % Target Start
    plot3(Actual_state_history(1, sel_Rtx), ...
          Actual_state_history(1, sel_Rtz), ...
          Actual_state_history(1, sel_Rty), ...
          'or', 'LineWidth', 2, 'MarkerFaceColor', 'r')
    % Missile Start
    plot3(Actual_state_history(1, sel_Rpx), ...
          Actual_state_history(1, sel_Rpz), ...
          Actual_state_history(1, sel_Rpy), ...
          'ob', 'LineWidth', 2, 'MarkerFaceColor', 'b')

    % 4. Correct Labels and View
    xlabel('Downrange [m]', 'FontSize', 16)
    ylabel('Crossrange [m]', 'FontSize', 16) % Visual Y is now Crossrange
    zlabel('Altitude [m]', 'FontSize', 16)   % Visual Z is now Altitude
    
    legend('Missile', 'Target', 'Target Start', 'Missile Start', 'Location', 'best')
    
    set(gca, 'fontsize', 16);
    set(gcf, 'color', 'w');
    grid on
    axis equal
    view(3) % Standard 3D view

    figure(2)
    plot(time(1:miss_index), aM(1:miss_index)./9.81, 'LineWidth', 2)
    xlabel('Time [s]', 'FontSize', 16)
    ylabel('Accelaration [g]', 'FontSize', 16)
    set(gca, 'fontsize', 16, 'ylim', [0 60]);
    set(gcf, 'color', 'w');
    grid on

    figure(3)
    plot(time, VC_hist, 'LineWidth', 2);
    xlabel('Time [s]', 'FontSize', 16)
    ylabel('Closing velocity [m/s]', 'FontSize', 16)
    set(gca, 'fontsize', 16);
    set(gcf, 'color', 'w');
    grid on

    figure(4)
    plot(time, vecnorm(LOS_rate, 2, 2)*180/pi, 'b', 'LineWidth', 2);
    xlabel('Time [s]', 'FontSize', 16)
    ylabel('LOS rate [deg/s]', 'FontSize', 16)
    set(gca, 'fontsize', 16);
    set(gcf, 'color', 'w');
    grid on

    figure(5)
    % Plot Magnitude of RTM
    semilogy(time, vecnorm(RTM, 2, 2), 'b', 'LineWidth', 2);
    xlabel('Time [s]', 'FontSize', 16)
    ylabel('RTM [m]', 'FontSize', 16)
    set(gca, 'fontsize', 16);
    set(gcf, 'color', 'w');
    grid on

end

if ~(strcmp(engagement_plots, 'on') || strcmp(engagement_plots, 'off'))
        disp('Not acceptable value. Enter ''on'' or ''off'' for engagement_plots');
end


% Animation
if strcmp(full_scale_engagement_viz, 'on')

    % 1. Setup Figure Window
    figure(6)
    clf;
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 800]);
    set(gcf, 'color', 'w');

    % 2. Calculate STATIC Limits (Global)
    % Swapping Y and Z so Altitude is Vertical (Z-axis in plot)
    all_X = [Actual_state_history(:, sel_Rpx); Actual_state_history(:, sel_Rtx)];
    all_Y = [Actual_state_history(:, sel_Rpz); Actual_state_history(:, sel_Rtz)]; % Visual Y is Crossrange
    all_Z = [Actual_state_history(:, sel_Rpy); Actual_state_history(:, sel_Rty)]; % Visual Z is Altitude

    margin = 0.05;
    xlim_fixed = [min(all_X), max(all_X)] + [-1, 1] * (range(all_X)*margin + 100);
    ylim_fixed = [min(all_Y), max(all_Y)] + [-1, 1] * (range(all_Y)*margin + 100);
    zlim_fixed = [min(all_Z), max(all_Z)] + [-1, 1] * (range(all_Z)*margin + 100);

    % 3. Apply Initial Settings
    title([PN_type ' ProNav, ' num2str(HE_pitch*180/pi) ' Heading error, N =' num2str(N)], 'FontSize', 16)
    xlabel('Downrange [m]', 'FontSize', 16)
    ylabel('Crossrange [m]', 'FontSize', 16)
    zlabel('Altitude [m]', 'FontSize', 16)
    
    grid on; axis equal; view(3);
    set(gca, 'fontsize', 16, 'xlim', xlim_fixed, 'ylim', ylim_fixed, 'zlim', zlim_fixed);
    hold on

    % 4. Plot Start Points
    plot3(Actual_state_history(1,sel_Rpx), Actual_state_history(1,sel_Rpz), Actual_state_history(1,sel_Rpy), ...
        'ob', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    plot3(Actual_state_history(1,sel_Rtx), Actual_state_history(1,sel_Rtz), Actual_state_history(1,sel_Rty), ...
        'or', 'LineWidth', 2, 'MarkerFaceColor', 'r');

    % 5. Animation Loop
    k=1;
    % FIX: Limit loop to available data size
    max_ind = min(miss_index, size(Actual_state_history, 1));

    for i = 1 : dt_index : max_ind
        
        if i > dt_index
            prev = i - dt_index;
            
            % Plot Segments (Swapping Y/Z)
            % Missile (Blue)
            plot3([Actual_state_history(i,sel_Rpx), Actual_state_history(prev,sel_Rpx)], ...
                [Actual_state_history(i,sel_Rpz), Actual_state_history(prev,sel_Rpz)], ...
                [Actual_state_history(i,sel_Rpy), Actual_state_history(prev,sel_Rpy)], ...
                'b-', 'LineWidth', 2);
            
            % Target (Red)
            plot3([Actual_state_history(i,sel_Rtx), Actual_state_history(prev,sel_Rtx)], ...
                [Actual_state_history(i,sel_Rtz), Actual_state_history(prev,sel_Rtz)], ...
                [Actual_state_history(i,sel_Rty), Actual_state_history(prev,sel_Rty)], ...
                'r-', 'LineWidth', 2);
        end

        drawnow; 
        % pause(0.01); % Uncomment if it's too fast

        if strcmp(vid_file, 'on')
            frame = getframe(gcf);
            % Force resize current frame to match the first frame exactly
            if k > 1
                frame.cdata = imresize(frame.cdata, size(F1(1).cdata, 1:2));
            end
            F1(k) = frame;
            k = k+1;
        end

    end
    
    if(strcmp(vid_file, 'on'))
        video = VideoWriter(['Media\Full_engagement_3D_' PN_type]);
        video.FrameRate = 30;
        open(video)
        writeVideo(video, F1)
        close(video);
    end
   
end

if ~(strcmp(full_scale_engagement_viz, 'on') || strcmp(full_scale_engagement_viz, 'off'))
        disp('Not acceptable value. Enter ''on'' or ''off'' for full_scale_engagement_viz');
end

if strcmp(zoomed_engagement_viz, 'on')

    k=1;
    title('');

    for i = 1 : dt_index : size(Actual_state_history, 1)

        figure(7)
        
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0, 0.5, 1])
        
        % Velocity Vector (3D)
        VMx = Actual_state_history(i, sel_Vpx);
        VMy = Actual_state_history(i, sel_Vpy);
        VMz = Actual_state_history(i, sel_Vpz);

        ph1 = quiver3(Actual_state_history(i, sel_Rpx), Actual_state_history(i, sel_Rpy), Actual_state_history(i, sel_Rpz), ...
            VMx, VMy, VMz, 'b', 'linewidth', 2, 'AutoScaleFactor', 1);
        hold on
        grid on
        axis equal
        view(3)

        if i >= 2
            % Trajectory tails (3D)
            plot3([Actual_state_history(i,sel_Rtx), Actual_state_history(i-dt_index,sel_Rtx)], ...
                [Actual_state_history(i,sel_Rty), Actual_state_history(i-dt_index,sel_Rty)], ...
                [Actual_state_history(i,sel_Rtz), Actual_state_history(i-dt_index,sel_Rtz)], ...
                'r-', 'LineWidth', 2)
            plot3([Actual_state_history(i,sel_Rpx), Actual_state_history(i-dt_index,sel_Rpx)], ...
                [Actual_state_history(i,sel_Rpy), Actual_state_history(i-dt_index,sel_Rpy)], ...
                [Actual_state_history(i,sel_Rpz), Actual_state_history(i-dt_index,sel_Rpz)], ...
                'b-', 'LineWidth', 2)
        end

        % Acceleration Vector Logic (3D Cross Product)
        % We calculate the direction using Cross(Omega, V)
        
        curr_Omega = LOS_rate(i, :);
        
        if strcmp(PN_type, 'Pure')
            curr_V = [VMx, VMy, VMz];
        elseif strcmp(PN_type, 'True')
            % Draw LOS line
            ph2 = plot3([Actual_state_history(i,sel_Rpx), Actual_state_history(i,sel_Rpx)+RTM(i,1)], ...
                        [Actual_state_history(i,sel_Rpy), Actual_state_history(i,sel_Rpy)+RTM(i,2)], ...
                        [Actual_state_history(i,sel_Rpz), Actual_state_history(i,sel_Rpz)+RTM(i,3)], ...
                'k--', 'LineWidth', 2);
            curr_V = -VTM(i,:); % Closing velocity vector
        else
            disp('PN_type must be ''True'' or ''Pure''');
        end
        
        % Calculate Acceleration Vector (Direction & Magnitude)
        acc_vec = cross(curr_Omega, curr_V) * N;
        
        % Scale acceleration for visualization
        if norm(acc_vec) > 0
             acc_vec = acc_vec / norm(acc_vec) * 1000; % Scale for visibility
        end

        ph3 = quiver3(Actual_state_history(i, sel_Rpx), Actual_state_history(i, sel_Rpy), Actual_state_history(i, sel_Rpz), ...
            acc_vec(1), acc_vec(2), acc_vec(3), ...
            'k', 'linewidth', 2);

        ph4 = plot3(Actual_state_history(i, sel_Rpx), Actual_state_history(i, sel_Rpy), Actual_state_history(i, sel_Rpz), ...
            'b.', 'LineWidth', 2, 'MarkerSize', 20);
        
        % Dynamic limits
        cx = Actual_state_history(i, sel_Rpx);
        cy = Actual_state_history(i, sel_Rpy);
        cz = Actual_state_history(i, sel_Rpz);

        set(gca, 'xlim', [cx-4e3, cx+4e3], ...
                 'ylim', [cy-4e3, cy+4e3], ...
                 'zlim', [cz-4e3, cz+4e3]);
                 
        set(gcf, 'color', 'w');
        title(['Engagement Visualization - ' PN_type ' ProNav'], 'FontSize', 14);

        pause(0.1);

        if strcmp(vid_file, 'on')
            F2(k) = getframe(gcf);
            k = k + 1;
        end

        if strcmp(PN_type, 'True')
            delete(ph2)
        end

        delete(ph1)
        delete(ph3)
        delete(ph4)

    end

    if(strcmp(vid_file, 'on'))
        video2 = VideoWriter(['Media\Zoomed_engagement_3D_' PN_type]);
        video2.FrameRate = 30;
        open(video2)
        writeVideo(video2, F2)
        close(video2);
    end

end

if ~(strcmp(zoomed_engagement_viz, 'on') || strcmp(zoomed_engagement_viz, 'off'))
        disp('Not acceptable value. Enter ''on'' or ''off'' for zoomed_engagement_viz');
end