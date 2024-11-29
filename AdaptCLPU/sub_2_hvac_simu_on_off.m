%% Coded by Rongxing Hu in May 2024, guliusha@outlook.com, NCSU
%% Coded and tested on Matlab 2019
%% set the HVAC status based on the scheduled status of each group
%% simulate the house HVAC load

function [Time_hvac_simu_plot_end] = sub_2_hvac_simu_on_off(step_time, d, otg_sh, hrzn_h,...
                                    input_fd, input_fd_pv, input_fd_load, input_fd_hvac, input_fd_ckt,...
                                    output_fd, output_fd_step,...
                                    k_res, k_pref, k_penalty_PVcurt, mark_ctr, mark_clpu,mark_case, mark_case_loop, mark_case_title_loop,...
                                    N_grp, N_house, M_house_idx_my_grp_ph,...
                                    pref_time, hvac_para, Tout_data_hour, ...
                                    EMSresult_each);

close all
time_calc_all_start_sub = tic;

%% set interval/ step
step_time = step_time;        % min; interval
step_time_hour = step_time/60; % in hour
N_step = hrzn_h*round(60/step_time); % steps, 

%% ==== select data length, day, hour
day1  = d ;
hour1 = otg_sh ;
step_start = ((day1-1)*24+hour1-1)*(60/step_time)+1;
step_end   = step_start + N_step - 1;

%% set import folder and output folder
input_fd = input_fd;
input_fd_pv = input_fd_pv; % load data
input_fd_load = input_fd_load; % load data
input_fd_hvac = input_fd_hvac; % 
input_fd_ckt = input_fd_ckt;    % circuit data

output_fd = output_fd;            % for EMS results 
output_fd_step = output_fd_step;  % for each EMS rolling results

output_fd_hvac_simu = strcat(output_fd,'\','Output_hvac_simu','\');  % for each EMS rolling results
if ~exist(output_fd_hvac_simu,'dir')
    mkdir(output_fd_hvac_simu)
end

%% set reserve / comfort parameters
k_res = k_res;                %  reserve factor, (1+k_res)*total_load<=total_gen
fprintf(['\n',' k_res = : ',num2str(k_res),'\n',])

k_pref     = k_pref;            %1.5  %1.1  Customer preferred served time period
fprintf(['\n',' k_pref = : ',num2str(k_pref),'\n',])

k_penalty_PVcurt = k_penalty_PVcurt;      % penalty factor in objective function for PV curtailment
fprintf(['\n',' penalty_PVcurt = : ',num2str(k_penalty_PVcurt),'\n',])

mark_ctr = mark_ctr;
mark_clpu = mark_clpu;

decay_delta_default = 0; % [0,1,2,3] must be integer;
fprintf(['\n',' decay_delta_default = : ',num2str(decay_delta_default),' (min)','\n',])

mark_case = mark_case; % from EMS
mark_case_simu = strcat(mark_case,'_Dc',num2str(decay_delta_default) );
mark_case_title_simu = replace(mark_case_simu,'_','-');

%% import circuit data
N_house = N_house;
N_grp = N_grp; %

M_house_idx_my_grp_ph = M_house_idx_my_grp_ph; % [grp_ph, house_idx] to mat; my house index list in per phase per group

%% import HVAC data that simulated using generated random parameters
hvac_para = hvac_para;
Phvac_rate_house = hvac_para.P_rate/1000; % in kW
Phvac_rate_house = Phvac_rate_house(1,1:N_house);
Phvac_rate_house = Phvac_rate_house';

%% data within  The selected otg period
pref_time = pref_time;
Tout_data_hour = Tout_data_hour;

%% import EMS-hvac results
EMS_result = EMSresult_each; % num2str(s),
sol_Pserve_hvac = EMS_result.sol_Pserve_grp_hvac;    % this is the group-ph load
sol_Uhouse = EMS_result.sol_Userve_hs;

sol_Uhouse_minu = kron(sol_Uhouse,ones(1,step_time));

%% get the HVAC status
hvac_staus = sol_Uhouse_minu; % in minu
hvac_simu_prepare.hvac_staus = hvac_staus;
 
save(strcat(output_fd_hvac_simu,  'hvac_simu_prepare','_',mark_case_simu,'.mat'), 'hvac_simu_prepare')

%% ======================================= run hvac simulation
Ca = hvac_para.Ca';
Cm = hvac_para.Cm';
R1 = hvac_para.R1';
R2 = hvac_para.R2';
Tset_all = hvac_para.Tset'; % in degree C
Phvac_rate_house = hvac_para.P_rate'; % in in W
Ca = Ca(1:N_house,:);
Cm = Cm(1:N_house,:);
R1 = R1(1:N_house,:);
R2 = R2(1:N_house,:);
Tset_all = Tset_all(1:N_house,:);
Phvac_rate_house = Phvac_rate_house(1:N_house,:);

TDB_base = hvac_para.TDB; % 1, in C degree

hour_add = 15; % add hours for initialize the simulation
Tout_hour_simu = [ones(hour_add,1)*Tout_data_hour(1);Tout_data_hour];
N_hour_simu = length(Tout_hour_simu);
tt = 1:N_hour_simu+1;
tx = 1:1/60:N_hour_simu+1;
Tout_hour_simu_original = [Tout_hour_simu;Tout_hour_simu(end)];
Tout_minu_simu = interp1(tt,Tout_hour_simu_original,tx)';
N_timestep = length(Tout_minu_simu);

Tout = Tout_minu_simu;
Timestep = 1;

house_pow_test = zeros(N_timestep,N_house);
house_status_test = zeros(N_timestep,N_house);
house_Troom_test = zeros(N_timestep,N_house); 

hvac_staus_simu = [zeros(N_house,hour_add*60),hvac_staus,hvac_staus(:,end)];

record_hvac_ini = zeros(N_house,4); %[Tset,x_ini,stat_ini]

for i = 1:N_house
    
    Tset_pick = Tset_all(i);
    x_ini0 = Tset_pick + TDB_base*(rand(1,1)-0.5); 
    x_ini  = [x_ini0';x_ini0'];  % [air temperature inside the house; mass temperature inside the house]
    stat_ini = randi([0,1],1,1); % initial status
    
    record_hvac_ini(i,:) = [Tset_pick,x_ini',stat_ini];

    A = [-(1/(R2(i)*Ca(i))+1/(R1(i)*Ca(i)))  1/(R2(i)*Ca(i)); 1/(R2(i)*Cm(i))  -1/(R2(i)*Cm(i))];
    B = [1/(R1(i)*Ca(i)) 1/Ca(i) ;0 0];
    C = [1 0; 0 1];
    D = [0; 0];
    P_house = Phvac_rate_house(i);
    
    house_ac_status = hvac_staus_simu(i,:);

    %% set initial values
    x0    = x_ini(:,1); % x_ini(:,i) = [x_ini0;x_ini0];
    stat0 = stat_ini(1);
    Tset = Tset_pick(1)*ones(N_timestep,1);
    TDB  = TDB_base*ones(N_timestep,1); % same deadband

    stat = zeros(N_timestep, 1);
    Q    =  zeros(N_timestep, 1);
    Troom_record = zeros(N_timestep,1);
    Tmass_record = zeros(N_timestep,1);
    
    if P_house <= 1000 || Tset_pick > 26 || Tset_pick < 16  % unreasonable HVAC parameter -> no HVAC
        fprintf(['**unreasonable HVAC parameter for house: ', num2str(i),'\n',])
    end

    %% --- pre-event, normal condition 
    N_step_normal = N_timestep;
    Tout_sub = Tout(1:N_timestep);

    [Troom_sub, Tmass_sub, stat_sub, Q_sub] = TCL_HVAC_simu_on_off(A, B, C, D, P_house, x0, stat0, Tset, TDB, Tout_sub, N_step_normal, Timestep, house_ac_status); % summer cooling mode

    house_pow_test(:,i) = Q_sub/1000; % in kW [N_step, N_house]
    house_status_test(:,i) = stat_sub;
    house_Troom_test(:,i) = Troom_sub;
    Troom_record(1:N_timestep,1) = Troom_sub;
    Tmass_record(1:N_timestep,1) = Tmass_sub;
%     % -- updated the inital temp and status
%     x0 = [Troom_sub(end);Tmass_sub(end)];
%     stat0 = stat_sub(end);

end 

%% save result % only save those within temperature profiles
% for houses
house_pow = house_pow_test(hour_add*60+1:N_hour_simu*60,:); % without the test period considering initial status
house_status = house_status_test(hour_add*60+1:N_hour_simu*60,:);
house_Troom = house_Troom_test(hour_add*60+1:N_hour_simu*60,:);

delta_tempC_sum_house = sum(house_Troom)/60 - N_step*step_time*Tset_all'/60 ; % compared to the Tset, in C*hour

Phvac_simu_house.house_pow  = house_pow;
Phvac_simu_house.house_status  = house_status;
Phvac_simu_house.house_Troom  = house_Troom;
Phvac_simu_house.temp_hour  = Tout_data_hour;
Phvac_simu_house.temp_minu  = Tout_minu_simu;
Phvac_simu_house.delta_tempC_sum_house = delta_tempC_sum_house;
Phvac_simu_house.record_hvac_ini = record_hvac_ini;

save(strcat(output_fd_hvac_simu,  'Phvac_simu_house_1m','_',mark_case_simu,'.mat'), 'Phvac_simu_house')

%% ========== get 60-min simulation data 
%--result record----------------- house
hvac_house_pow_60_ori = zeros(N_house,N_step); % average
N_hour = N_step;

for n = 1:N_house
    Phouse = house_pow(:,n);
    Qp_60 = reshape(Phouse,60,N_hour*60/60);
    P_house_60 = mean(Qp_60);        % in kW 
    hvac_house_pow_60_ori(n,:) = P_house_60;
end

%% == save result for each resol
Phvac_simu_house_60m.hvac_house_pow  = hvac_house_pow_60_ori';
Phvac_simu_house_60m.temp_hour  = Tout_data_hour;
Phvac_simu_house_60m.temp_minu  = Tout_minu_simu(hour_add*60+1:N_hour_simu*60);

save(strcat(output_fd_hvac_simu,  'Phvac_simu_house_60m','_',mark_case_simu,'.mat'), 'Phvac_simu_house_60m')

%% plot result
plot_hvac_simu_RT

Time_hvac_simu_plot_end = toc(time_calc_all_start_sub);
fprintf(['\n','######## Finished HVAC simu & plots In: ',num2str(Time_hvac_simu_plot_end),' sec.', '\n',])