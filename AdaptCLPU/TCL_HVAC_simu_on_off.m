%% Coded by Rongxing Hu in May 2024, guliusha@outlook.com, NCSU
%% Coded and tested on Matlab 2019
%% sub-func @ House HVAC load simulation @ 1-min 
% summer always in cooling mode

function [Troom_sub, Tmass_sub, stat_sub, Q_sub, T_on,T_off] = TCL_HVAC_simu_on_off(A, B, C, D,P_house, x0, stat0, Tset, TDB, Tout_sub, N_step_normal, Timestep, house_ac_status)

N_timestep= N_step_normal;
Tout = [Tout_sub; Tout_sub(end)];           % +1
Tset = [Tset; Tset(end,1)]; % +1
TDB  = [TDB; TDB(end,1)];    % +1
mark_cycle_control_step = 0*[house_ac_status, house_ac_status(end)];    % +1
house_ac_status = [house_ac_status, house_ac_status(end)];

stat = zeros(N_timestep+1, 1);
Q    =  zeros(N_timestep+1, 1);
Troom_record = zeros(N_timestep+1,1);
Tmass_record = zeros(N_timestep+1,1);

T_on = zeros(N_timestep+1,1);
T_off = zeros(N_timestep+1,1);

%-- set event status, binary
U_house_ac = house_ac_status;

% first step
stat(1) = stat0*U_house_ac(1);
if stat(1) == 1
    Q(1) = P_house;    % *****
    T_on(1,1) = 1;
    T_off(1,1) = 0;
else
    Q(1) = 0;
    T_on(1,1) = 0;
    T_off(1,1) = 1;
end

cool_mode = 1;

if cool_mode == 0 %% heating mode
    u0 = [Tout(1); Q(1)];
else
    u0 = [Tout(1); -1*Q(1)];
end

% the other steps
for j = 2:1:N_timestep+1 % 
    DT = ( A*x0 + B*u0 )*Timestep*60; % minute to second, linearization
    x1 = x0 + DT;
    
    if mark_cycle_control_step(j) == 0   % not in cycle-control mode
        % cooling mode
        if x1(1)>(Tset(j) + TDB(j)/2) &&  stat(j-1) == 0
            stat(j) = 1*U_house_ac(j);
        elseif x1(1)< (Tset(j) - TDB(j)/2) &&  stat(j-1) == 1
            stat(j) = 0;
        else
            stat(j) = stat(j-1)*U_house_ac(j);
        end
        if stat(j) == 1
            Q(j) = P_house*U_house_ac(j);% *Timestep; % ***** W*s
            T_on(j,1) = T_on(j-1,1) + 1;
            T_off(j,1) = 0;
        else
            Q(j) = 0;
            T_on(j,1) = 0;
            T_off(j,1) = T_off(j-1,1) + 1;
        end
        x0 = x1;
        u0 = [Tout(j); -1*Q(j)];
        Troom_record(j-1) = x1(1);
        Tmass_record(j-1) = x1(2);
    else
        % cooling mode
        if x1(1)>(Tset(j) + TDB(j)/2) &&  stat(j-1) == 0
            stat(j) = 1*U_house_ac(j);
        elseif x1(1)< (Tset(j) - TDB(j)/2) &&  stat(j-1) == 1
            stat(j) = 0;
        else
            stat(j) = stat(j-1)*U_house_ac(j);
        end
        if stat(j) == 1
            Q(j) = P_house*U_house_ac(j);% *Timestep; % ***** W*s
            T_on(j,1) = T_on(j-1,1) + 1;
            T_off(j,1) = 0;
        else
            Q(j) = 0;
            T_on(j,1) = 0;
            T_off(j,1) = T_off(j-1,1) + 1;
        end
        x0 = x1;
        u0 = [Tout(j); -1*Q(j)];
        Troom_record(j-1) = x1(1);
        Tmass_record(j-1) = x1(2);

    end
                        
end
DT = A*x0  + B*u0; %the last second
x1 = x0 + DT;
Troom_record(N_timestep) = x1(1);
Tmass_record(N_timestep) = x1(2);

Troom_sub = Troom_record;
Tmass_sub = Tmass_record;
stat_sub = stat;
Q_sub = Q;

Troom_sub = Troom_record(1:end-1,1);
Tmass_sub = Tmass_record(1:end-1,1);

stat_sub = stat(1:end-1,1);
Q_sub = Q(1:end-1,1);
T_on = T_on(1:end-1,1);
T_off = T_off(1:end-1,1);


