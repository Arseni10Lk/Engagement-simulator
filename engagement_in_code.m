clc
clear all
close all

% dbstop if naninf


% Visualisation data
engagement_plots = 'on';
full_scale_engagement_viz = 'off';
zoomed_engagement_viz = 'off';

vid_file = 'off';

% initial conditions
PN_type = 'True';
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
for j = 1 : length_time-1

    s1 = nlinpronav_sim(Actual_state, N, at, PN_type);
    s2 = nlinpronav_sim(Actual_state + t_step*s1/2, N, at, PN_type);
    s3 = nlinpronav_sim(Actual_state + t_step*s2/2, N, at, PN_type);
    s4 = nlinpronav_sim(Actual_state + t_step*s3, N, at, PN_type);

    Actual_state = Actual_state + t_step*(s1+2*s2+2*s3+s4)/6;

    Con_vs_t(:, j + 1) = Actual_state;
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

% Visualization
dt_index = .1/t_step;
miss_index = find(min(abs(RTM)) == abs(RTM));

if strcmp(engagement_plots, 'on')

    figure(1)
    plot(Actual_state(1:miss_index ,sel_Rpx), Actual_state(1:miss_index ,sel_Rpy), 'LineWidth', 2)
    hold on 
    plot(Actual_state(1:miss_index,sel_Rtx), Actual_state(1:miss_index,sel_Rty), 'r--', 'LineWidth', 2)
    plot(Actual_state(1,sel_Rtx), Actual_state(1,sel_Rty), 'or', 'LineWidth', 2)
    plot(Actual_state(1,sel_Rpx), Actual_state(1,sel_Rpy), 'ob', 'LineWidth', 2)
    xlabel('Downrange [m]', 'FontSize', 16)
    ylabel('Altitude [m]', 'FontSize', 16)
    ylim([0 15000])
    set(gca, 'fontsize', 16);
    set(gcf, 'color', 'w');
    grid on


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
    set(gca, 'fontsize', 16, 'ylim', [4000 5000]);
    set(gcf, 'color', 'w');
    grid on

    figure(4)
    plot(time, lambda_dot*180/pi, 'b', 'LineWidth', 2);
    xlabel('Time [s]', 'FontSize', 16)
    ylabel('LOS rate [m/s]', 'FontSize', 16)
    set(gca, 'fontsize', 16, 'ylim', [-5 5]);
    set(gcf, 'color', 'w');
    grid on

    figure(5)
    semilogy(time, RTM, 'b', 'LineWidth', 2);
    xlabel('Time [s]', 'FontSize', 16)
    ylabel('RTM [m]', 'FontSize', 16)
    set(gca, 'fontsize', 16, 'xlim', [7 9]);
    set(gcf, 'color', 'w');
    grid on

end

if ~(strcmp(engagement_plots, 'on') || strcmp(engagement_plots, 'off'))

        disp('Not acceptable value. Enter ''on'' or ''off'' for engagement_plots');

end


% Animation
if strcmp(full_scale_engagement_viz, 'on')

    k=1;
    for i = 1 : dt_index : (miss_index + dt_index)
        if i == 1
            
            figure(6)
            plot(Actual_state(1,sel_Rpx), Actual_state(1,sel_Rpy), 'ob', 'LineWidth', 1, 'MarkerFaceColor', 'b')
            hold on
            plot(Actual_state(1,sel_Rtx), Actual_state(1,sel_Rty), 'or', 'LineWidth', 1, 'MarkerFaceColor', 'r')
            title([PN_type ' ProNav, ' num2str(HE_pitch*180/pi) ' Heading error, N =' num2str(N)], 'FontSize', 16)
            xlabel('Downrange [m]', 'FontSize', 16)
            ylabel('Altitude [m]', 'FontSize', 16)
            set(gca, 'fontsize', 16, 'xlim', [0 40000], 'ylim', [6000 12000], ...
                'position', [0.1220 0.1381 0.8388 0.7119])
            set(gcf, 'color', 'w', 'position', [10, 80, 1212, 298]);
            grid on
            axis equal

        end

        if i >= 2
            
            figure(6)

            plot([Actual_state(i,sel_Rpx), Actual_state(i-dt_index,sel_Rpx)], ...
                [Actual_state(i,sel_Rpy), Actual_state(i-dt_index,sel_Rpy)], ...
                'b-', 'LineWidth', 2)
            hold on
            plot([Actual_state(i,sel_Rtx), Actual_state(i-dt_index,sel_Rtx)], ...
                [Actual_state(i,sel_Rty), Actual_state(i-dt_index,sel_Rty)], ...
                'r-', 'LineWidth', 2)
            set(gca, 'fontsize', 16, 'xlim', [0 40000], 'ylim', [6000 12000], ...
                'position', [0.1220 0.1381 0.8388 0.7119]);
        end

        pause(0.1)

        if strcmp(vid_file, 'on')
            F1(k) = getframe(gcf);
            k = k+1;
        end

    end
    
    if(strcmp(vid_file, 'on'))
        video = VideoWriter(['Media\Full_engagement_20deg_HE_' PN_type]);
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

    for i = 1 : dt_index : (miss_index + dt_index)

        figure(7)
        
        VMx = Actual_state(i, sel_Vpx);
        VMy = Actual_state(i, sel_Vpy);

        ph1 = quiver(Actual_state(i, sel_Rpx), Actual_state(i, sel_Rpy), VMx, VMy, ...
            'b', 'linewidth', 2);
        hold on
        grid on

        if i >= 2
            
            plot([Actual_state(i,sel_Rtx), Actual_state(i-dt_index,sel_Rtx)], ...
                [Actual_state(i,sel_Rty), Actual_state(i-dt_index,sel_Rty)], ...
                'r-', 'LineWidth', 2)
            plot([Actual_state(i,sel_Rpx), Actual_state(i-dt_index,sel_Rpx)], ...
                [Actual_state(i,sel_Rpy), Actual_state(i-dt_index,sel_Rpy)], ...
                'b-', 'LineWidth', 2)

        end

        if strcmp(PN_type, 'Pure')

            Heading_pursuer = atan2(VMy, VMx);
            aMx = -aM(i) * sin(Heading_pursuer)*8;
            aMy = aM(i) * cos(Heading_pursuer)*8;

        elseif strcmp(PN_type, 'True')
            
            ph2 = plot([Actual_state(i,sel_Rpx), Actual_state(i,sel_Rpx)+RTMx(i)], ...
                [Actual_state(i,sel_Rpy), Actual_state(i,sel_Rpy)+RTMy(i)], ...
                'k--', 'LineWidth', 2);

            aMx = -aM(i) * sin(lambda_pitch(i))*8;
            aMy =  aM(i) * cos(lambda_pitch(i))*8;

        else

            disp('PN_type must be ''True'' or ''Pure''');

        end

        ph3 = quiver(Actual_state(i, sel_Rpx), Actual_state(i, sel_Rpy), aMx, aMy, ...
            'k', 'linewidth', 2);

        ph4 = plot(Actual_state(i, sel_Rpx), Actual_state(i, sel_Rpy), ...
            'b.', 'LineWidth', 2, 'MarkerSize', 20);

        set(gca, 'xlim', [Actual_state(i, sel_Rpx)-4e3, Actual_state(i, sel_Rpx)+4e3], ...
            'ylim', [Actual_state(i, sel_Rpy)-4e3, Actual_state(i, sel_Rpy)+4e3]);
        set(gcf, 'color', 'w', 'position', [360 278 560 420]); % ?
        title(['Engagement Visualization - ' PN_type ' ProNav'], ...
            'FontSize', 14);

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
        video2 = VideoWriter(['Media\Zoomed_engagement_20deg_HE_' PN_type]);
        video2.FrameRate = 30;
        open(video2)
        writeVideo(video2, F2)
        close(video2);
    end

end

if ~(strcmp(zoomed_engagement_viz, 'on') || strcmp(zoomed_engagement_viz, 'off'))

        disp('Not acceptable value. Enter ''on'' or ''off'' for zoomed_engagement_viz');

end

