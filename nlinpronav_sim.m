function dy = nlinpronav_sim(t, y, HE, Np, aT, VM, ProNav_type)

sel_beta = 1;
sel_RT1 = 2;
sel_RT2 = 3;
sel_RM1 = 4;
sel_RM2 = 5;
sel_VT1 = 6;
sel_VT2 = 7;
sel_VM1 = 8;
sel_VM2 = 9;

dy = [0;0;0;0;0;0;0;0;0];

% collecting data from input
VT = sqrt(y(sel_VT1)^2+y(sel_VT2)^2);

RTM1 = y(sel_RT1) - y(sel_RM1);
RTM2 = y(sel_RT2) - y(sel_RM2);
RTM = sqrt(RTM1^2+RTM2^2);

VTM1 = y(sel_VT1) - y(sel_VM1);
VTM2 = y(sel_VT2) - y(sel_VM2);

lambda = atan2(RTM2, RTM1);
lambda_dot = (RTM1*VTM2 - RTM2*VTM1)/RTM^2;

VC = -(RTM1*VTM1 + RTM2*VTM2)/RTM;

% calculation for return values

dy(1) = aT/VT;
dy(2) = y(sel_VT1);
dy(3) = y(sel_VT2);
dy(4) = y(sel_VM1);
dy(5) = y(sel_VM2);
dy(6) = aT*sin(y(sel_beta));
dy(7) = aT*cos(y(sel_beta));

if strcmp(ProNav_type, 'True')
    nc = Np*VC*lambda_dot;
    dy(8) = -nc*sin(lambda);
    dy(9) = nc*cos(lambda);
elseif strcmp(ProNav_type, 'Pure')
    Vp_horizon_angle = atan2(y(sel_VM2), y(sel_VM1));
    nc = Np*VM*lambda_dot;
    dy(8) = -nc*sin(Vp_horizon_angle);
    dy(9) = nc*cos(Vp_horizon_angle);
else
    disp("Wrong type");
    return
end

end