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

%% modelling and solving

%%

time_finish_all = toc(time_Calc_all_start);
fprintf(['\n','######## Finished In: ',num2str(time_finish_all),' sec.', '\n',]) 