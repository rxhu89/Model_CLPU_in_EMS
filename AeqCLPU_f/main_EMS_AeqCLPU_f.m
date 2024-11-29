%% Coded by Rongxing Hu in May 2024, guliusha@outlook.com, NCSU
%% Coded and tested on Matlab 2019
%% EMS for testing CLPU models @ AeqCLPU
%% using yalmip for modeling EMS problems, using Gurobi 11 for solving the MILP EMS problem
%% load: nonHVAC & normal HVAC load and CLPU part

clc; clear all; close all

time_Calc_all_start = tic;
start_clock = datestr(datetime('now'));
fprintf(['\n', '***** Running the code at: ', start_clock, ' *****\n'])

%% set interval/ step, data resol, days, outages
step_time = 60;                % in min; interval
step_time_hour = step_time/60; % in hour

%% Set Outage / opti horizon --------
hrzn_h = 48;                         % in hour
N_step = hrzn_h*round(60/step_time); % total steps, 
fprintf([' stg_1 horizon (h): ',num2str(hrzn_h),'\n',])

otg_sh = 1;  % the start hour for each run

N_days_all_ori = 24;
if otg_sh == 1
    N_days_start_all = N_days_all_ori - (hrzn_h/24) + 1; % start-days
else
    N_days_start_all = N_days_all_ori - (hrzn_h/24); % start-days
end

%% set input_folder and output_folder
dictory_parts = strsplit(pwd, '\');
dictory_input_parent = strjoin(dictory_parts(1:end-1), '\');

input_fd = strcat(dictory_input_parent,'\Input_data\');
input_fd_pv = strcat(input_fd,'\Input_pv','\');   % pv
input_fd_load = strcat(input_fd,'\Input_load','\'); % Non-HVAC load data
input_fd_hvac = strcat(input_fd,'\Input_hvac\');  % HVACs
input_fd_ckt = strcat(input_fd,'\Input_ckt\');    % circuit data

%% set reserve / comfort parameters
k_res = 0.15;                %  reserve factor, (1+k_res)*total_load <= total_gen
fprintf(['\n',' k_res = : ',num2str(k_res),'\n',])

k_pref     = 1.5;            %  Customer preferred served time period
fprintf(['\n',' k_pref = : ',num2str(k_pref),'\n',])

k_penalty_PVcurt = 0.0001;   % penalty factor in objective function for PV curtailment
fprintf(['\n',' penalty_PVcurt = : ',num2str(k_penalty_PVcurt),'\n',])

Dpeak = 2; % if CLPU, it denotes the max. peak duration
fprintf(['\n',' Dpeak(Max.) = : ',num2str(Dpeak),' hours','\n',])

mark_ctr = strcat('grp');    % mark load-serving control
mark_clpu = strcat(' CLPU Peak (vary) @ k_clpu @ multiple steps @ hard constrained');
fprintf(['\n',' Case setting = : ',mark_ctr,'\n',])
fprintf(['\n',' CLPU setting = : ',mark_clpu,'\n',])

% mark the preferred periods @ just set for all days-----
pref_time0 = ones(1, 24*(60/step_time) );
priod_t1   = [8,9];   % -8 -9 
priod_t2   = [19,20]; % -19 -20

pref_time0(1,(priod_t1(1)-1)*(60/step_time)+1:priod_t1(2)*(60/step_time)) = ...
                      k_pref*pref_time0(1,(priod_t1(1)-1)*(60/step_time)+1:priod_t1(2)*(60/step_time) ); % the preferred period 7-9 18-20
pref_time0(1,(priod_t2(1)-1)*(60/step_time)+1:priod_t2(2)*(60/step_time)) = ...
                      k_pref*pref_time0(1,(priod_t2(1)-1)*(60/step_time)+1:priod_t2(2)*(60/step_time) );
pref_time0 = repmat(pref_time0,1,N_days_all_ori); %  to a month

mark_case = strcat(mark_ctr,'_','kpref_',num2str(k_pref),'_','Dp_m', num2str(Dpeak));
mark_title = replace(mark_case,'_','-');

%% set the components
P_PV_rate  = 4500/3 ;         % kW, per ph

P_bat_rate = 3000/3;          % kW, per ph 
E_bat_rate = P_bat_rate*3*2;  % kWh, 3-phase, max(PV_data)/2*2;
k_dc   =  0.95;               % discharge and charge efficency
% operational settings
soc_min = 0.2;
soc_max = 0.9;

%% import circuit data
M_house_idx_my_grp_ph = importdata(strcat(input_fd_ckt,'\','map_house_idx_grp_ph','.mat'));

N_house = sum(M_house_idx_my_grp_ph(:,:) > 0,'all');
vec_Nhouse_per_grp_ph = sum(M_house_idx_my_grp_ph(:,:) > 0,2); % [N_grp_ph,1]

%% import topol candidate mapping data
map_grp_topol = readtable( strcat(input_fd_ckt,'IEEE123_topology_candidates_reduced.xlsx'),'sheet','group_topol','ReadRowNames',true,'ReadVariableNames',true);
map_sw_topol = readtable( strcat(input_fd_ckt,'IEEE123_topology_candidates_reduced.xlsx'),'sheet','sw_topol','ReadRowNames',true,'ReadVariableNames',true);
map_grp_topol = map_grp_topol{:,:}; % [N_grp,N_topol]
map_sw_topol = map_sw_topol{:,:}; % [N_sw,N_topol]

N_topol = size(map_sw_topol,2);
N_sw = size(map_sw_topol,1);
N_grp = size(map_grp_topol,1);

%% import PV dataset
Ppv_data_ori = importdata( strcat(input_fd_pv,'\','PV_data','.mat') ); % in kW
% PV_rate_raw = Ppv_data_ori.Ppv_rate;
Ppv_data0 = Ppv_data_ori.Ppv_data; % in 60-min resol

%% import load data (Non-HVAC)
Pgrp_ph_Nonhvac_data0 =  importdata( strcat(input_fd_load,'\','Pgrp_ph_Nonhvac_data','.mat') ); % in 60-min resol; [N_grp_ph, N_step] % Aug data
Pgrp_ph_Nonhvac_data0(isnan(Pgrp_ph_Nonhvac_data0)) = 0;

%% import temperature profiles
Tout_data_hour0 = importdata(strcat(input_fd_hvac,'\','Tout_data_hour_C.mat'));% in Celsius; [N_step,1]

%% import HVAC data that simulated using generated random parameters
hvac_para = importdata( strcat(input_fd_hvac,'\','hvac_para_random_2000_35_3000_Sduty1','.mat') ); 
Phvac_rate_house = hvac_para.P_rate/1000; % in kW
Phvac_rate_house = Phvac_rate_house(1,1:N_house);
Phvac_rate_house = Phvac_rate_house';

%% HVAC normal load data
Pgrp_ph_hvac_norm_data0 = importdata( strcat(input_fd_hvac,'\','Pgrp_ph_hvac_norm_data','.mat') ); % [N_grp_ph, N_step]

%% to group data format
Phvac_rate_grp_ph = zeros(N_grp*3, 1);
for i = 1:N_grp
    for j = 1:3
        temp_house_idx_ph_grp = M_house_idx_my_grp_ph((i-1)*3+j,:);
        temp_house_idx_ph_grp(temp_house_idx_ph_grp==0) = []; % delete zeros
        
        Phvac_rate_grp_ph((i-1)*3+j,1) = sum(Phvac_rate_house(temp_house_idx_ph_grp,1));
    end
end
Pgrp_ph_norm_data0 = Pgrp_ph_Nonhvac_data0 + Pgrp_ph_hvac_norm_data0;

%% set k_clpu @ changing how to calculate the Pclpu & Kclpu to get NoCLPU model and FixCLPU model(also needs to adjust the Dpeak)
Pclpu_grp_ph_data0 = Phvac_rate_grp_ph - Pgrp_ph_hvac_norm_data0; %% here uses this approach to calculate Pclpu/Kclpu !!

k_clpu_grp_data0 = zeros(N_grp, length(Tout_data_hour0)); %[N_grp, N_step]; the clpu factor @ fixed  Phvac = Phvac_norm + Pclpu
for i = 1:N_grp
    temp_Pclpu = Pclpu_grp_ph_data0((i-1)*3+1:i*3,:);
    temp_Phvac_rate = Phvac_rate_grp_ph((i-1)*3+1:i*3,:);
    temp_kclpu_grp = sum(temp_Pclpu,1)/sum(temp_Phvac_rate);
    
    for j = 1:length(Tout_data_hour0)
        temp_temperature = round(Tout_data_hour0(j)); % only get the parameter table with 1 C interval
        if temp_temperature < 22 % if increasing/decreasing, case by case?
            temp_kclpu_grp(1,j) = 0;
        end
        if temp_temperature > 40
            fprintf(['!!! temaperature > 40 C'])
            temp_kclpu_grp(1,j) =  0; % HVACs are working at fully power at high temperature
        end
    end
    k_clpu_grp_data0(i,:) = temp_kclpu_grp;
end

%% set Dpeak for each hours @ updated
Tout_data_ave_proc0 = Tout_data_hour0; % if temperature is decreasing, get the ave_temp of the last 2 steps
for i = 1:length(Tout_data_hour0)
    if i == 1
        Tout_data_ave_proc0(i) = Tout_data_hour0(i);
    elseif i == 2
        if Tout_data_hour0(i-1) > Tout_data_hour0(i)
            Tout_data_ave_proc0(i) = Tout_data_hour0(i-1);
        end
    else
        if Tout_data_hour0(i-2) > Tout_data_hour0(i) && Tout_data_hour0(i-1) > Tout_data_hour0(i)
            Tout_data_ave_proc0(i) = mean(Tout_data_hour0(i-2:i-1));
        end
    end
end

Dpeak_grp_data0 = ones(N_grp, length(Tout_data_hour0));       
for i = 1:length(Tout_data_ave_proc0)
    temp_temperature = round(Tout_data_ave_proc0(i)); %
    if temp_temperature >= 32 && temp_temperature < 37 % 32
        Dpeak_grp_data0(:,i) = Dpeak * ones(N_grp,1);
    elseif  temp_temperature >= 37
        Dpeak_grp_data0(:,i) = Dpeak * ones(N_grp,1);
    end
end           

%% set some inital values
soc_0 = 0.9;
Userve_grp_0 = zeros(N_grp,1); % initial servced status of groups
dre_grp_0 = zeros(N_grp,1); % initial remaining peak duration

M_big = N_step; % only for CLPU consumption calculation

for d = 1 : 1 % N_days_start_all
%% rolling placehold
    close all;
    
    time_calc_step_start = tic; % get the calculation time of the step

    %% ==== select data length, day, hour
    day1  = d ;
    hour1 = otg_sh ;

    step_start = ((day1-1)*24+hour1-1)*(60/step_time)+1;
    step_end   = step_start + N_step - 1;

    fprintf(['\n****',' day_1 = : ',num2str(day1),' Hour_start = ',num2str(hour1),' duration (h) = ',num2str(hrzn_h),]);
    %--
    temp_mark_period = strcat('d',num2str(day1),'_','sh',num2str(hour1));
    mark_case_loop = strcat(mark_case,'_',temp_mark_period);
    mark_case_title_loop = replace(mark_case_loop,'_','-');
    
    fprintf(['\n','!!!!!!!!! @ Starting Case @@ ',mark_case_loop, '\n',])
    
    %% == set output_folder
    output_data_fd_name = strcat( 'hz',num2str(hrzn_h),'_','d',num2str(day1),'_sh',num2str(hour1) );
    output_fd = strcat(pwd,'\Output_Data\',output_data_fd_name); % for EMS results 
    output_fd_step = strcat(output_fd,'\','each_step');  % 
    if ~exist(output_fd,'dir')
        mkdir(output_fd)
    end
    if ~exist(output_fd_step,'dir')
        mkdir(output_fd_step)
    end
    
    %% == data within the selected optimization horizon / otg period
    Ppv_data = Ppv_data0(step_start:step_end);

    Pgrp_ph_Nonhvac_data = Pgrp_ph_Nonhvac_data0(:,step_start:step_end);
    Pgrp_ph_hvac_norm_data = Pgrp_ph_hvac_norm_data0(:,step_start:step_end);

    pref_time = pref_time0(step_start:step_end);
    Tout_data_hour = Tout_data_hour0(step_start:step_end);
    %-- CLPU
    k_clpu_grp = k_clpu_grp_data0(:,step_start:step_end);
    Dpeak_grp_data = Dpeak_grp_data0(:,step_start:step_end);
    
    %% initial values
    soc_ini = soc_0;
    Userve_grp_ini = Userve_grp_0;
    dre_grp_ini = dre_grp_0;

    %% ***************** modelling **************************************
    %% == define all variables =
    yalmip('clear')  % https://github.com/yalmip/YALMIP/discussions/956; yalmip doesn't clear the previous same-name variables in a loop like matlab

    x_Ppv   = sdpvar(1,N_step,'full');  % PV power
    x_Qpv   = sdpvar(1,N_step,'full');  % PV power
    x_Ppv_curt  = sdpvar(1,N_step,'full');  % PV curtailment
    
    x_Ebat     = sdpvar(1,N_step,'full'); % BESS energy
    x_Pc    = sdpvar(3,N_step,'full');  % BESS charging power
    x_Pd    = sdpvar(3,N_step,'full');  % BESS discharging power
    x_Pbat  = sdpvar(3,N_step,'full');  % BESS power
    x_Ud    = binvar(3,N_step,'full');  % BESS discharging status
    x_Qbat    = sdpvar(3,N_step,'full');% BESS reactive power

    x_Pserve_grp  = sdpvar(N_grp*3,N_step,'full'); % served load P; per group per phase
    x_Qserve_grp  = sdpvar(N_grp*3,N_step,'full'); % served load Q; per group per phase
    x_Pserve_pcc  = sdpvar(3,N_step,'full'); % served total load P on phase; PCC
    x_Qserve_pcc  = sdpvar(3,N_step,'full'); % served total load P on phase; PCC 
    x_Pserve_fd = sdpvar(1,N_step,'full');   % total served load in each interval, PCC total

    x_Pserve_grp_hvac = sdpvar(N_grp*3,N_step,'full'); % served load P; per group per phase

    x_Utopol  = binvar(N_topol,N_step,'full');   % topology status
    x_Userve_grp  = binvar(N_grp,N_step,'full'); % group serving status
    x_Usw     = binvar(N_sw,N_step,'full');      % switch status

    x_Uclpu_grp  = binvar(N_grp,N_step,'full');   % clpu event status
    x_Pclpu_grp  = sdpvar(N_grp*3,N_step,'full'); % clpu consumption

    x_dpeak_grp = sdpvar(N_grp,N_step,'full'); % peak dutration; obtained at the first restoration step
    x_Ure_grp  = binvar(N_grp,N_step,'full');  % auxiliary variable for calculating the remaing CLPU duration
    x_dre_grp  = sdpvar(N_grp,N_step,'full');  % remainging CLPU (peak) duration

    %% ==== contstraints =====
    constr = [];

    %% PV
    constr = [constr, ( 0 <= x_Ppv <= Ppv_data ):'PV_limit_P'];
    constr = [constr, ( 0 <= x_Qpv <= P_PV_rate ):'PV_limit_Q'];
    constr = [constr, ( x_Ppv + x_Qpv <= sqrt(2)*P_PV_rate ):'PV_limit_OP'];
    constr = [constr, ( x_Ppv_curt == Ppv_data -x_Ppv ):'PV_curt_P'];

    %% battery
    constr = [constr, ( soc_min*E_bat_rate <= x_Ebat <= soc_max*E_bat_rate ):'Bat_E_lb'];
    %-- 
    constr = [constr, ( x_Ebat(1) == soc_ini*E_bat_rate + sum( (k_dc*x_Pc(:,1) - x_Pd(:,1)/k_dc), 1 )*step_time_hour ):'Bat_E_balan1' ];
    if N_step >1
        constr = [constr, ( x_Ebat(2:end) == x_Ebat(1:end-1) + sum( (k_dc*x_Pc(:,2:end) - x_Pd(:,2:end)/k_dc),1 )*step_time_hour ):'Bat_E_balan2' ]; 
    end

    constr = [constr, ( 0 <= x_Pc <= P_bat_rate*(1-x_Ud) ):'Bat_ch' ];
    constr = [constr, ( 0 <= x_Pd <= P_bat_rate*x_Ud ):'Bat_dch' ];
    constr = [constr, ( x_Pbat == x_Pd - x_Pc ):'Bat_pow' ];   % PCC interface

    constr = [constr, ( 0 <= x_Qbat <= P_bat_rate ):'Bat_Q1']; % may have oversized inverter?
    constr = [constr, ( x_Pbat + x_Qbat <= sqrt(2)*P_bat_rate ):'Bat_Q2']; 
    constr = [constr, ( -1*x_Pbat + x_Qbat <= sqrt(2)*P_bat_rate ):'Bat_Q3'];

    %% load serving
    for g = 1:N_grp
        for p = 1:3
            constr = [constr, ( x_Pserve_grp((g-1)*3+p,:) == ...
                      x_Userve_grp(g,:) .* Pgrp_ph_Nonhvac_data((g-1)*3+p,:) +  x_Pserve_grp_hvac((g-1)*3+p,:) ):['serve_Pgrp' num2str(g) '_'  num2str(p)]];
        end
    end
    constr = [constr, ( x_Pserve_fd == sum(x_Pserve_grp,1) ):'serve_Pfd' ];
    for p = 1:3
        constr = [constr, ( x_Pserve_pcc(p,:) == sum(x_Pserve_grp(p:3:end, :) ) ):['serve_Pph_' num2str(p)] ];
    end
    constr = [constr, ( x_Qserve_grp == 0.5*x_Pserve_grp ):'serve_Qhs' ];
    constr = [constr, ( x_Qserve_pcc == 0.5*x_Pserve_pcc ):['serve_Qph'] ];

    %% HVAC @ CLPU
    constr_hvac = [];
    for g = 1:N_grp
        for p = 1:3
            constr_hvac = [constr_hvac, ( x_Pserve_grp_hvac((g-1)*3+p,:) == ...
                           x_Userve_grp(g,:) .* Pgrp_ph_hvac_norm_data((g-1)*3+p,:) + x_Pclpu_grp((g-1)*3+p,:) ):['serve_Pgrp_hvac' num2str(g) '_'  num2str(p)]];
            constr_hvac = [constr_hvac, ( x_Pclpu_grp((g-1)*3+p,:) == ...
                           k_clpu_grp(g,:) .* x_Uclpu_grp(g,:) .* Phvac_rate_grp_ph((g-1)*3+p,:) ):['serve_Pgrp_clpu' num2str(g) '_'  num2str(p)]];
        end
    end

    constr_hvac = [constr_hvac, ( 0 <= x_dpeak_grp <= Dpeak_grp_data .* x_Userve_grp ):'dpeak_1' ];
    constr_hvac = [constr_hvac, ( x_dpeak_grp(:,1) >= Dpeak_grp_data(:,1) .* (x_Userve_grp(:,1)-Userve_grp_ini) ):'dpeak_2' ];
    constr_hvac = [constr_hvac, ( x_dpeak_grp(:,2:end) >= Dpeak_grp_data(:,2:end) .* (x_Userve_grp(:,2:end)-x_Userve_grp(:,1:end-1)) ):'dpeak_3' ];
    constr_hvac = [constr_hvac, ( x_dpeak_grp(:,1) <= Dpeak_grp_data(:,1) .* (1-Userve_grp_ini) ):'dpeak_4' ];
    constr_hvac = [constr_hvac, ( x_dpeak_grp(:,2:end) <= Dpeak_grp_data(:,2:end) .* (1-x_Userve_grp(:,1:end-1)) ):'dpeak_5' ];

    constr_hvac = [constr_hvac, ( x_dre_grp >= 0 ):'dre_1' ];
    constr_hvac = [constr_hvac, ( x_dre_grp(:,1) >= dre_grp_ini + x_dpeak_grp(:,1) - Userve_grp_ini - Dpeak_grp_data(:,1) .* (1-x_Userve_grp(:,1)) ):'dre_2' ];
    constr_hvac = [constr_hvac, ( x_dre_grp(:,2:end) >= x_dre_grp(:,1:end-1) + x_dpeak_grp(:,2:end) ...
                                                        - x_Userve_grp(:,1:end-1) - Dpeak_grp_data(:,2:end) .* (1-x_Userve_grp(:,2:end)) ):'dre_3' ];

    constr_hvac = [constr_hvac, ( x_dre_grp <= 0 + M_big*(1-x_Ure_grp) ):'dre_4' ];
    constr_hvac = [constr_hvac, ( x_dre_grp(:,1) <= dre_grp_ini + x_dpeak_grp(:,1) - Userve_grp_ini ...
                                                    - Dpeak_grp_data(:,1) .* (1-x_Userve_grp(:,1)) + M_big*x_Ure_grp(:,1) ):'dre_4' ];
    constr_hvac = [constr_hvac, ( x_dre_grp(:,2:end) <= x_dre_grp(:,1:end-1) + x_dpeak_grp(:,2:end) - x_Userve_grp(:,1:end-1) ...
                                                        - Dpeak_grp_data(:,2:end) .* (1-x_Userve_grp(:,2:end)) + M_big*x_Ure_grp(:,2:end) ):'dre_5' ];

    constr_hvac = [constr_hvac, ( x_Uclpu_grp >= x_dre_grp ./ Dpeak_grp_data ):'Uclpu_1' ]; %% note error when Any Dpeak == 0
    constr_hvac = [constr_hvac, ( x_Uclpu_grp <= x_dre_grp ):'Uclpu_2' ];

    %% balance
    for p = 1:3
        constr = [constr, ( x_Pserve_pcc(p,:) == x_Ppv + x_Pbat(p,:) ):['Balan_P_' num2str(p)] ];
        constr = [constr, ( x_Qserve_pcc(p,:) == x_Qpv + x_Qbat(p,:) ):['Balan_Q_' num2str(p)] ];
    end

    %% reserve
    constr_res = [];
    for p = 1:3
        constr_res = [constr_res, ( (1+k_res) * x_Pserve_pcc(p,:) <= Ppv_data + P_bat_rate ):['Reserve_P_' num2str(p)] ];
    end
    constr_res = [constr_res, ( k_res*x_Pserve_fd*step_time_hour <= 3*x_Ppv_curt*step_time_hour + (E_bat_rate-x_Ebat) ):['Reserve_E_'] ];

    %% topology
    constr_tp = [];
    for i = 1:N_step
        constr_tp = [constr_tp, ( x_Userve_grp(:,i) == map_grp_topol*x_Utopol(:,i) ):['Topol_Ug_' num2str(i)] ];
        constr_tp = [constr_tp, ( x_Usw(:,i) == map_sw_topol*x_Utopol(:,i) ):['Topol_Usw_' num2str(i)] ];    
    end
    constr_tp = [constr_tp, ( sum(x_Utopol,1) == ones(1,N_step) ):'Topol_Utopol' ];

    %% objective function
    f_Etotal_pref = sum(x_Pserve_fd.*pref_time) * step_time_hour;
    f_penalty_PVcurt = k_penalty_PVcurt*sum(x_Ppv_curt) * step_time_hour * 3 ; % all phases

    objective = -1*f_Etotal_pref + f_penalty_PVcurt;

    %% =======solve =============================================
    options = sdpsettings('solver','gurobi','debug',1,'verbose',1);  % verbose: display,0: no; 1 simple-display; 2 full-display
    options.savesolverinput = 1;   % 
    options.savesolveroutput = 1;  % 
    options.gurobi.TuneTimeLimit = 0; % for Gurobi error 10008: Unable to set parameter TuneTimeLimit to value -1 (minimum is 0), check https://groups.google.com/g/yalmip/c/2VuGmw-0oOI
    
    solver_time_max_set = 60*10;   % in sec
    options.gurobi.timelimit = solver_time_max_set; 

    mip_gap = 0.0001; % default 0.0001 @ relative tolerance on the gap between the best integer objective 
    options.gurobi.mipgap = mip_gap; % works 

    fprintf(['\n','******** Beginning solve ','\n',])

    sol = optimize([constr,constr_hvac,constr_res,constr_tp],objective,options);

    fval = value(objective);
    fprintf(['\n','******** fval: ',num2str(fval),'(kWh)\n',])

    %% Analyze error flags
    mark_get_sol = 0;
    if sol.problem == 0 || sol.problem == 3 % 0 Successfully solved; 1 Infeasible problem; 2 Unbounded objective function; 3 Maximum iterations exceeded
        redu_sd = 0;
        fprintf(['******** fval: ',num2str(fval),'\n',])
        fprintf(['   ^^^^ ',sol.info,'\n']) 
        fprintf(['^^^^ sol_MIPgap: ',num2str(sol.solveroutput.result.mipgap),'\n']) % gurobi can do this; CPLEX cannot
        mark_get_sol = 1;
    else  
        disp('Hmm, something went wrong!');    
        fprintf(['\n',' ^^^^ ',sol.info,'\n'])      
        error('!!!! Failure !!!! ')
    end

    %% get the results
    sol_info = sol;

    sol_Ppv = value(x_Ppv);  
    sol_Qpv = value(x_Qpv);  
    sol_Ppv_curt = value(x_Ppv_curt);  
    sol_Ebat = value(x_Ebat);
    sol_Pc   = value(x_Pc);  
    sol_Pd   = value(x_Pd);  
    sol_Pbat = value(x_Pbat);  
    sol_Ud   = value(x_Ud);  
    sol_Qbat = value(x_Qbat);

    sol_Pserve_grp = value(x_Pserve_grp); 
    sol_Qserve_grp = value(x_Qserve_grp); 
    sol_Pserve_pcc = value(x_Pserve_pcc); 
    sol_Qserve_pcc = value(x_Qserve_pcc);  
    sol_Pserve_fd = value(x_Pserve_fd); 

    sol_Pserve_grp_hvac = value(x_Pserve_grp_hvac); 

    sol_Utopol = value(x_Utopol); 
    sol_Userve_grp = value(x_Userve_grp);
    sol_Usw = value(x_Usw); 

    sol_Uclpu_grp = value(x_Uclpu_grp); 
    sol_Pclpu_grp = value(x_Pclpu_grp); 

    sol_dpeak_grp  = value(x_dpeak_grp);
    sol_Ure_grp  = value(x_Ure_grp);
    sol_dre_grp  = value(x_dre_grp);

    %-
    sol_Ud = round(sol_Ud);
    sol_Utopol = round(sol_Utopol);
    sol_Userve_grp = round(sol_Userve_grp);
    sol_Usw = round(sol_Usw);
    sol_Uclpu_grp = round(sol_Uclpu_grp);
    sol_Ure_grp  = round(sol_Ure_grp);

    %% processed 
    temp_sol_Userve_grp = kron(sol_Userve_grp,ones(3,1));
    sol_Pserve_grp_Nonhvac = temp_sol_Userve_grp .* Pgrp_ph_Nonhvac_data;
    sol_Pserve_grp_hvac_norm = temp_sol_Userve_grp .* Pgrp_ph_hvac_norm_data;

    %% get house status
    sol_Userve_hs = zeros(N_house,N_step);
    for g = 1:N_grp
        temp_house_idx_in_grp = M_house_idx_my_grp_ph( (g-1)*3+1:g*3,:);
        temp_house_idx_in_grp = unique(temp_house_idx_in_grp);
        temp_house_idx_in_grp(temp_house_idx_in_grp==0) = [];

        sol_Userve_hs(temp_house_idx_in_grp,:) = repmat(sol_Userve_grp(g,:),length(temp_house_idx_in_grp),1);
    end

    save_sol_mipgap = sol.solveroutput.result.mipgap; % gurobi can do this; CPLEX cannot

    time_calc_step = toc(time_calc_step_start);
    fprintf(['---- Solved in : ',num2str(time_calc_step),' sec.','\n'])

    %% save results of each rolling
    sub_1_save_ems_sol_grp_clpu_XpdF

    %% plot EMS results
    plot_EMS_grp_DA_clpu
    
    %% run HVAC simulations
    [Time_hvac_simu_plot_end] = sub_2_hvac_simu_on_off(step_time, d, otg_sh, hrzn_h,...
                                input_fd, input_fd_pv, input_fd_load, input_fd_hvac, input_fd_ckt,...
                                output_fd, output_fd_step,...
                                k_res, k_pref, k_penalty_PVcurt, Dpeak, mark_ctr, mark_clpu, mark_case, mark_case_loop, mark_case_title_loop,...
                                N_grp, N_house, M_house_idx_my_grp_ph,...
                                pref_time, hvac_para, Tout_data_hour, ...
                                EMSresult_each);

end

time_finish_all = toc(time_Calc_all_start);
fprintf(['\n','######## Finished In: ',num2str(time_finish_all),' sec.', '\n',]) 