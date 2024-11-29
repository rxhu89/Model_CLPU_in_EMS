%% Coded by Rongxing Hu in May 2024, guliusha@outlook.com, NCSU
%% Coded and tested on Matlab 2019
%% save EMS results 

EMSresult_each.k_res = k_res;
EMSresult_each.k_pref = k_pref;
EMSresult_each.k_penalty_PVcurt = k_penalty_PVcurt;
% result_each.Dpeak = Dpeak;

EMSresult_each.data_period = [hrzn_h, day1, otg_sh]; %
EMSresult_each.pref_time = pref_time;

EMSresult_each.sys_component.P_PV_rate = P_PV_rate;
EMSresult_each.sys_component.P_bat_rate = P_bat_rate;
EMSresult_each.sys_component.E_bat_rate = E_bat_rate;
EMSresult_each.sys_component.k_dc = k_dc;
EMSresult_each.sys_component.soc_min = soc_min;
EMSresult_each.sys_component.soc_max = soc_max;

EMSresult_each.opti_ini.soc_0 = soc_0;
EMSresult_each.opti_ini.Userve_grp_0 = Userve_grp_0;
EMSresult_each.opti_ini.clpu_Dpeak_satu_0 = clpu_Dpeak_satu_0;
EMSresult_each.opti_ini.clpu_dp_grp_0 = clpu_dp_grp_0;
EMSresult_each.opti_ini.clpu_dpeak_grp_0 = clpu_dpeak_grp_0;
EMSresult_each.opti_ini.clpu_dre_grp_0 = clpu_dre_grp_0;
EMSresult_each.opti_ini.clpu_khvac_grp_0 = clpu_khvac_grp_0;

EMSresult_each.Ppv_data = Ppv_data;
EMSresult_each.Pgrp_ph_Nonhvac_data = Pgrp_ph_Nonhvac_data;
EMSresult_each.Pgrp_ph_hvac_norm_data = Pgrp_ph_hvac_norm_data;
EMSresult_each.Q_over_p = 0.5;
%-- clpu
EMSresult_each.Tout_data_hour = Tout_data_hour;
EMSresult_each.Phvac_rate_grp_ph = Phvac_rate_grp_ph;
EMSresult_each.clpu_kpeak_dura = clpu_kpeak_dura;
EMSresult_each.clpu_kdecay_dura = clpu_kdecay_dura;
EMSresult_each.clpu_Dpeak_satu = clpu_Dpeak_satu;
EMSresult_each.clpu_knorm = clpu_knorm;
EMSresult_each.clpu_kpk = clpu_kpk;

EMSresult_each.sol_info = sol_info;
EMSresult_each.sol_Ppv = sol_Ppv;
EMSresult_each.sol_Qpv = sol_Qpv;
EMSresult_each.sol_Ppv_curt = sol_Ppv_curt;
EMSresult_each.sol_Ebat = sol_Ebat;
EMSresult_each.sol_Pc = sol_Pc;
EMSresult_each.sol_Pd = sol_Pd;
EMSresult_each.sol_Pbat = sol_Pbat;
EMSresult_each.sol_Ud = sol_Ud;
EMSresult_each.sol_Qbat = sol_Qbat;

EMSresult_each.sol_Pserve_grp = sol_Pserve_grp;
EMSresult_each.sol_Qserve_grp = sol_Qserve_grp;
EMSresult_each.sol_Pserve_pcc = sol_Pserve_pcc;
EMSresult_each.sol_Qserve_pcc = sol_Qserve_pcc;
EMSresult_each.sol_Pserve_fd = sol_Pserve_fd;

EMSresult_each.sol_Pserve_grp_hvac  = sol_Pserve_grp_hvac;

EMSresult_each.sol_Utopol = sol_Utopol;
EMSresult_each.sol_Userve_grp = sol_Userve_grp;
EMSresult_each.sol_Usw    = sol_Usw;

% result_each.sol_Uclpu_grp = sol_Uclpu_grp; 
% result_each.sol_Pclpu_grp = sol_Pclpu_grp; 
%--
% result_each.sol_dpeak_grp  = sol_dpeak_grp;
% result_each.sol_Ure_grp  = sol_Ure_grp;
% result_each.sol_dre_grp  = sol_dre_grp;
EMSresult_each.sol_dp_grp     = sol_dp_grp;
EMSresult_each.sol_dpeak_grp  = sol_dpeak_grp;
EMSresult_each.sol_dre_grp    = sol_dre_grp;
EMSresult_each.sol_Usatu_grp  = sol_Usatu_grp;
EMSresult_each.sol_Udecay_grp = sol_Udecay_grp;

EMSresult_each.sol_khvac_grp = sol_khvac_grp;
EMSresult_each.sol_kclpu_grp = sol_kclpu_grp;
EMSresult_each.sol_Pclpu_grp = sol_Pclpu_grp; 

% processed
EMSresult_each.sol_Pserve_grp_Nonhvac = sol_Pserve_grp_Nonhvac;    
EMSresult_each.sol_Pserve_grp_hvac_norm = sol_Pserve_grp_hvac_norm;

EMSresult_each.sol_Userve_hs = sol_Userve_hs;
EMSresult_each.sol_kclpu_hs = sol_kclpu_hs;

EMSresult_each.sol_time = time_calc_step;
EMSresult_each.save_sol_mipgap = save_sol_mipgap;

save(strcat(output_fd_step,'\','EMSresult_each','_',mark_case,'.mat'), 'EMSresult_each') % num2str(s),