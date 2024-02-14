clc
clear all
close all

% dbstop if naninf


% Visualisation data
engagement_plots = 'off';
full_scale_engagement_viz = 'on';
zoomed_engagement_viz = 'on';

vid_file = 'on';


% initial conditions
PN_type = 'Pure';
at = 0;
HE = -20*pi/180; % rad
N = 3;
tf = 12;
t_step = 0.01;

beta = 0; % rad
Rtx = 40000; % m
Rty = 10000; % m
Rpx = 0; % m
Rpy = 10000; % m

VP = 4000; % m/s
Vt = 1000; % m/s

% calculation of dependent inital conditions
Vtx = -Vt*cos(beta);
Vty = Vt*sin(beta);

RTMx = Rtx - Rpx;
RTMy = Rty - Rpy;

lambda = atan2(RTMy, RTMx); % rad

Lead_angle = asin(Vt*sin(beta + lambda)/VP); % rad

Vpx = VP*cos(Lead_angle + lambda + HE); % m/s 
Vpy = VP*sin(Lead_angle + lambda + HE); % m/s


% Vector with Initial conditions

y0 = [beta, Rtx, Rty, Rpx, Rpy, Vtx, Vty, Vpx, Vpy]';

% application of ProNav laws
t = 0:t_step:tf-t_step; lt = length(t);

Con_vs_t = zeros(length(y0), lt);
Con_vs_t(:, 1) =  y0; 

% Integration with RK4
Actual_state = y0;
Con_vs_t = Actual_state;
for j = 1 : lt-1
    s1 = nlinpronav_sim(t(j), Actual_state, HE, N, at, VP, PN_type);
    s2 = nlinpronav_sim(t(j)+t_step/2, Actual_state + t_step*s1/2, HE, N, at, VP, PN_type);
    s3 = nlinpronav_sim(t(j)+t_step/2, Actual_state + t_step*s2/2, HE, N, at, VP, PN_type);
    s4 = nlinpronav_sim(t(j) + t_step, Actual_state + t_step*s3, HE, N, at, VP, PN_type);
    Actual_state = Actual_state + t_step*(s1+2*s2+2*s3+s4)/6;
    Con_vs_t = [Con_vs_t, Actual_state];
end
Actual_state = Con_vs_t';

% Processing simulation results

sel_beta = 1;
sel_Rtx = 2;
sel_Rty = 3;
sel_Rpx  = 4;
sel_Rpy = 5;
sel_Vtx = 6;
sel_Vty = 7;
sel_Vpx = 8;
sel_Vpy = 9;

VP = sqrt(Actual_state(:, sel_Vpx).^2+Actual_state(:,sel_Vpy).^2);
Vt = sqrt(Actual_state(:, sel_Vtx).^2+Actual_state(:,sel_Vty).^2);

RTMx = Actual_state(:, sel_Rtx)-Actual_state(:, sel_Rpx);
RTMy = Actual_state(:, sel_Rty)-Actual_state(:, sel_Rpy);
VTM1 = Actual_state(:, sel_Vtx)-Actual_state(:, sel_Vpx);
VTM2 = Actual_state(:, sel_Vty)-Actual_state(:, sel_Vpy);

RTM = sqrt(RTMy.^2+RTMx.^2);

lambda = atan2(RTMy, RTMx);
lambda_dot = (RTMx.*VTM2 - RTMy.*VTM1)./RTM.^2;

VC = -(RTMx.*VTM1 + RTMy.*VTM2)./RTM;

if strcmp(PN_type, 'True')
    aM = N*VC.*lambda_dot;
elseif strcmp(PN_type, 'Pure')
    aM = N*VP.*lambda_dot;
else
    disp("Wrong type");
    return
end

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
    plot(t(1:miss_index), aM(1:miss_index)./32.2, 'LineWidth', 2)
    xlabel('Time [s]', 'FontSize', 16)
    ylabel('Accelaration [g]', 'FontSize', 16)
    set(gca, 'fontsize', 16, 'ylim', [0 20]);
    set(gcf, 'color', 'w');
    grid on

    figure(3)
    plot(t, VC, 'LineWidth', 2);
    xlabel('Time [s]', 'FontSize', 16)
    ylabel('Closing velocity [m/s]', 'FontSize', 16)
    set(gca, 'fontsize', 16, 'ylim', [4000 5000]);
    set(gcf, 'color', 'w');
    grid on

    figure(4)
    plot(t, lambda_dot*180/pi, 'b', 'LineWidth', 2);
    xlabel('Time [s]', 'FontSize', 16)
    ylabel('LOS rate [m/s]', 'FontSize', 16)
    set(gca, 'fontsize', 16, 'ylim', [-5 5]);
    set(gcf, 'color', 'w');
    grid on

    figure(5)
    semilogy(t, RTM, 'b', 'LineWidth', 2);
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
            title([PN_type ' ProNav, ' num2str(HE*180/pi) ' Heading error, N =' num2str(N)], 'FontSize', 16)
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

            aMx = -aM(i) * sin(lambda(i))*8;
            aMy =  aM(i) * cos(lambda(i))*8;

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

