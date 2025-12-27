function dy = nlinpronav_sim(y, Np, aT, ProNav_type) % y should be a column

sel_RT = 1:3;
sel_RM = 4:6;
sel_VT = 7:9;
sel_VM = 10:12;

dy = zeros(12, 1);

% collecting data from input

RTM = y(sel_RT) - y(sel_RM);

VTM = y(sel_VT) - y(sel_VM);

VC = -dot(RTM,VTM)/norm(RTM); % closing speed

% defining vectors for angle processing

if strcmp(ProNav_type, 'True')
    % For True ProNav, the reference is the Line-of-Sight
    ref_vec = RTM; 
    vel_ref = VC;
elseif strcmp(ProNav_type, 'Pure')
    % For Pure ProNav, the reference is Velocity
    ref_vec = y(sel_VM);
    vel_ref = norm(y(sel_VM));
else
    error('Wrong ProNav_type');
end

u_fwd = ref_vec / norm(ref_vec);

global_up = [0; 1; 0];

if norm(abs(u_fwd) - global_up) < 10e-6
    u_right = [0; 0; 1];
else 
    u_right = cross(u_fwd, global_up);
    u_right = u_right/norm(u_right);
end

u_local_up = cross(u_right, u_fwd);
u_local_up = u_local_up/norm(u_local_up);

R_mat = [u_right, u_local_up, u_fwd];

LOS_rate = cross(RTM, VTM)/norm(RTM)^2;
lambda_dot_pitch = dot(LOS_rate, u_right);
lambda_dot_yaw = dot(LOS_rate, u_local_up);

% Standard ProNav: a = N * V * lambda_dot
    
ac_pitch = Np * vel_ref * lambda_dot_pitch;
ac_yaw   = Np * vel_ref * lambda_dot_yaw;
ac_local = [ac_yaw; ac_pitch; 0];

% Global coordinates

ac_global = R_mat * ac_local;

% calculation for return values

dy(1:3) = y(sel_VT);
dy(4:6) = y(sel_VM);
dy(7:9) = aT(:);
dy(10:12) = ac_global;

end