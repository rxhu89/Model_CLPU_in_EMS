%% Coded by Rongxing Hu in May 2024, guliusha@outlook.com, NCSU
%% Coded and tested on Matlab 2019
%% plot EMS results

close all

%--- set figure size---------------------
T_otg = N_step;
hrzn_d = hrzn_h / 24;

vec_time_label = [step_time_hour : step_time_hour : 24/step_time_hour];
vec_time_label = repmat(vec_time_label, 1, (hrzn_d+1));
if otg_sh == 1
    temp_vec = [0 : step_time_hour : otg_sh/step_time_hour];
    vec_time_label = [temp_vec(1:end-1), vec_time_label(1:hrzn_h/step_time_hour)];
else
    vec_time_label = vec_time_label((otg_sh-1)/step_time_hour : step_time_hour : ((otg_sh-1)+hrzn_h)/step_time_hour);
end

if hrzn_d <= 1 
    fig_size   = [0.3,0.5,0.3,0.4]; % one day data
    fig_size_s = [0.3,0.4,0.3,0.45]; % for served load box
    fig_size_b = [0.3,0.4,0.3,0.5]; % for served load box
    xTick_time = [0:4:T_otg]./1;
    xTick_label = vec_time_label(1+[0:4:T_otg]);
elseif hrzn_d >1 && hrzn_d <=2
    fig_size   = [0.3,0.5,0.4,0.4]; % for 2-3 day data; served node is different
    fig_size_s = [0.3,0.4,0.5,0.45]; % for served load box
    fig_size_b = [0.3,0.4,0.4,0.5]; % for served load box
    xTick_time = [0:4:T_otg]./1;
    xTick_label = vec_time_label(1+[0:4:T_otg]);
elseif hrzn_d >=3 && hrzn_d <=4
    fig_size   = [0.2,0.5,0.6,0.4]; % for 4-5 day data; served node is different
    fig_size_s = [0.2,0.4,0.7,0.45]; % for served load box
    fig_size_b = [0.2,0.4,0.6,0.5]; % for served load box
    xTick_time = [0:6:T_otg]./1;
    xTick_label = vec_time_label(1+[0:6:T_otg]);
else
    fig_size   = [0.2,0.5,0.7,0.4]; % for 6-7 day data; served node is different
    fig_size_s = [0.1,0.4,0.8,0.45]; % for served load box
    fig_size_b = [0.1,0.4,0.7,0.5]; % for served load box
    xTick_time = [0:6:T_otg]./1;
    xTick_label = vec_time_label(1+[0:6:T_otg]);
end

%% data process
mark_pref = pref_time;
mark_pref(mark_pref==1) = 0;
mark_pref(mark_pref>1) = 1;

Pgrp_ph_norm_data = Pgrp_ph_Nonhvac_data + Pgrp_ph_hvac_norm_data; % the total normal load;
feeder_Ppv_data = 3*Ppv_data;

feeder_Ptotal_norm_demand = sum(Pgrp_ph_norm_data,1);
feeder_P_Nonhvac_demand = sum(Pgrp_ph_Nonhvac_data,1);
feeder_Phvac_norm_demand = sum(Pgrp_ph_hvac_norm_data,1);

feeder_P_Nonhvac_demand_ph = [Pgrp_ph_Nonhvac_data(1:3:end,:); Pgrp_ph_Nonhvac_data(2:3:end,:); Pgrp_ph_Nonhvac_data(3:3:end,:)];
feeder_Phvac_norm_demand_ph = [Pgrp_ph_hvac_norm_data(1:3:end,:); Pgrp_ph_hvac_norm_data(2:3:end,:); Pgrp_ph_hvac_norm_data(3:3:end,:)];

feeder_Ptotal_norm_demand_ph = feeder_P_Nonhvac_demand_ph + feeder_Phvac_norm_demand_ph;

E_norm_demand = sum(sum(Pgrp_ph_norm_data))*step_time_hour;
E_Nonhvac_demand = sum(sum(Pgrp_ph_Nonhvac_data))*step_time_hour;
Ehvac_norm_demand = sum(sum(Pgrp_ph_hvac_norm_data))*step_time_hour;

P_all_served = sum(sol_Pserve_grp_Nonhvac + sol_Pserve_grp_hvac,1);
P_all_served_norm = sum(sol_Pserve_grp_Nonhvac + sol_Pserve_grp_hvac_norm,1);

E_served_demand = sum(P_all_served)*step_time_hour;
E_pref_served = sum(P_all_served.*mark_pref)*step_time_hour;

Ehvac_all_served = sum(sum(sol_Pserve_grp_hvac))*step_time_hour;
Ehvac_all_served_norm = sum(sum(sol_Pserve_grp_hvac_norm))*step_time_hour;
Ehvac_all_served_clpu = sum(sum(sol_Pclpu_grp))*step_time_hour; % only clpu addition

E_Nonhvac_served = sum(sum(sol_Pserve_grp_Nonhvac))*step_time_hour;

Epv_fore = 3*sum(sum(Ppv_data))*step_time_hour; % all 3-ph
Epv_used = 3*sum(sum(sol_Ppv))*step_time_hour; % all 3-ph
Epv_curt = Epv_fore - Epv_used;

RowNames_str = {'E_norm_demand','E_Nonhvac_demand','Ehvac_norm_demand',...
                'E_served_demand','E_pref_served',...
                'Ehvac_all_served','Ehvac_norm_served','Ehvac_clpu_served',...
                'E_Nonhvac_served','Epv_fore','Epv_used','Epv_curt',};
summary_array = [E_norm_demand,E_Nonhvac_demand,Ehvac_norm_demand,...
                 E_served_demand, E_pref_served,...
                 Ehvac_all_served,Ehvac_all_served_norm,Ehvac_all_served_clpu,...
                 E_Nonhvac_served, Epv_fore, Epv_used, Epv_curt]';
table_summary = array2table(summary_array,'RowNames',RowNames_str );
writetable(table_summary,strcat(output_fd,'\','EMS_table_summary','_',mark_case,'.xlsx'),'WriteRowNames',true,'WriteVariableNames',false) ;


%% -- plot inputs

fig_input = figure;
set(gcf,'unit','normalized','position',fig_size); % the size of the figure)
plot([1:T_otg]/1,sum(Pgrp_ph_norm_data,1),'k-','LineWidth',1)
bar([1:T_otg]/1,[feeder_P_Nonhvac_demand;feeder_Phvac_norm_demand],'stack') % ,1 ) %orange [0.85,0.33,0.10]
hold on
plot([1:T_otg]/1,feeder_Ppv_data,'-','Color',[0.39,0.83,0.07],'LineWidth',2)

ylabel('kW')
xlabel('Hour')
% ylim([0 1.2])
title([' Feeder demand and PV',' @ ', mark_title ])
legend('NonHVAC load','HVAC load (norm)','PV','Location','northwest')   % % ,'Location','northwest','Orientation','vertical') %'horizontal'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',1) set the axis & title & legend font
set(gca,'xTick',xTick_time) % set(gca,'xTick',[0:4:T_otg]./2)
set(gca,'xTickLabel',xTick_label)
% set(gcf,'Color',[1,1,1])
saveas(fig_input,strcat(output_fd, '\','input','_',mark_case,'.fig'));


%% --- plot the overall
fig_perf = figure;
set(gcf,'unit','normalized','position',fig_size); % the size of the figure)
plot([1:T_otg]/1,sum(Pgrp_ph_norm_data,1),'k-','LineWidth',1)
hold on
plot([1:T_otg]/1,feeder_Ppv_data,'-','LineWidth',1)
hold on
plot([1:T_otg]/1,sum(sol_Pserve_grp,1),'r:','LineWidth',1.5)
hold on
this_bar = bar([1:T_otg]/1,[sum(sol_Pserve_grp_Nonhvac,1);sum(sol_Pserve_grp_hvac,1);],'stack'); % ,1 ) %orange [0.85,0.33,0.10]
this_bar(1).FaceColor = [0.60,0.87,0.24];
this_bar(2).FaceColor = [0.18,0.62,0.91];
% this_bar(3).FaceColor = [0.65,0.65,0.65];
hold on
plot([1:T_otg]/1,sol_Ppv*3,'--','LineWidth',1)
hold on
plot([1:T_otg]/1,sum(sol_Pbat,1),'--','LineWidth',1)

ylabel('kW')
xlabel('Hour')
% ylim([0 1.2])
title([' Performance',' @ ', mark_title ])
legend('Feeder load(norm)','PV forecast','Served load','Served NonHVAC','Served hvac','Utilized PV','BESS','Location','northwest')   % % ,'Location','northwest','Orientation','vertical') %'horizontal'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',1) set the axis & title & legend font
set(gca,'xTick',xTick_time) %
set(gca,'xTickLabel',xTick_label)
% set(gcf,'Color',[1,1,1])
saveas(fig_perf,strcat(output_fd, '\','perform','_',mark_case,'.fig'));

%% --- plot the load;
fig_perf_load = figure;
set(gcf,'unit','normalized','position',fig_size); % the size of the figure)
plot([1:T_otg]/1,sum(Pgrp_ph_norm_data,1),'-','LineWidth',1,'color',[0.85,0.33,0.10])
hold on
plot([1:T_otg]/1,sum(Pgrp_ph_Nonhvac_data,1),'k-','LineWidth',1)
hold on
plot([1:T_otg]/1,sum(Pgrp_ph_hvac_norm_data,1),'-','LineWidth',1)
hold on
plot([1:T_otg]/1,sum(sol_Pserve_pcc,1),'r:','LineWidth',1.5)
hold on
this_bar = bar([1:T_otg]/1,[sum(sol_Pserve_grp_Nonhvac,1);sum(sol_Pserve_grp_hvac,1);],1.2);% ,'stack') % ,1 ) %orange [0.85,0.33,0.10]
this_bar(1).FaceColor = [0.60,0.87,0.24];
this_bar(2).FaceColor = [0.18,0.62,0.91];
% this_bar(3).FaceColor = [0.65,0.65,0.65];

ylabel('kW')
xlabel('Hour')
% ylim([0 1.2])
title([' Performance',' @ ', mark_title ])
legend('Feeder load(norm)','NonHVAC load','Normal HVAC','Served load','Served NonHVAC','Served HVAC','Served other','Location','northwest')   % % ,'Location','northwest','Orientation','vertical') %'horizontal'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',1) set the axis & title & legend font
set(gca,'xTick',xTick_time) % set(gca,'xTick',[0:4:T_otg]./2)
set(gca,'xTickLabel',xTick_label)
% set(gcf,'Color',[1,1,1])
saveas(fig_perf_load,strcat(output_fd, '\','perf_load','_',mark_case,'.fig'));

%% --- plot Battery
fig_bess = figure;
set(gcf,'unit','normalized','position',fig_size); % the size of the figure)

yyaxis left
plot([1:T_otg]/1,-1*sum(sol_Pc,1),'k-','LineWidth',1)
hold on
plot([1:T_otg]/1,sum(sol_Pd,1),'r-','LineWidth',1)
hold on
plot([1:T_otg]/1,sum(sol_Pbat,1),'g:','LineWidth',1.5)
hold on
plot([1:T_otg]/1,sum(sol_Qbat,1),'--','LineWidth',1)
ylabel('kW (dis+)')
xlabel('Hour')

yyaxis right
plot([1:T_otg]/1,sol_Ebat,'-','LineWidth',1.5)
ylim([0,E_bat_rate])
ylabel('kWh')

% ylim([0 1.2])
title(['Battery',' @ ', mark_title ])
legend('Pch','Pdis','Pbat','Qbat','Ebat','Location','northwest')   % % ,'Location','northwest','Orientation','vertical') %'horizontal'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',1) set the axis & title & legend font
set(gca,'xTick',xTick_time) % set(gca,'xTick',[0:4:T_otg]./2)
set(gca,'xTickLabel',xTick_label)
% set(gcf,'Color',[1,1,1])
saveas(fig_bess,strcat(output_fd, '\','bess','_',mark_case,'.fig'));

%% --- plot PCC
fig_pcc = figure;
set(gcf,'unit','normalized','position',fig_size); % the size of the figure)

hplot = plot([1:T_otg]/1,sol_Pserve_pcc,'LineWidth',1);
set(hplot, {'color'}, {[1.0,0.7,0.0];[0 1 0]; [1 0 0]})

ylabel('kW ')
xlabel('Hour')

% ylim([0 1.2])
title(['PCC',' @ ', mark_title ])
legend('a', 'b', 'c','Location','northwest')   % % ,'Location','northwest','Orientation','vertical') %'horizontal'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',1) set the axis & title & legend font
set(gca,'xTick',xTick_time) % set(gca,'xTick',[0:4:T_otg]./2)
set(gca,'xTickLabel',xTick_label)
% set(gcf,'Color',[1,1,1])
saveas(fig_pcc,strcat(output_fd, '\','pcc','_',mark_case,'.fig'));

%% -----plot served time of each group------------
%% other version: the best one
S1_fg_Ug = figure; 
set(gcf,'unit','normalized','position',[fig_size(1:3),0.5*fig_size(4)]); % the size of the figure
imagesc(sol_Userve_grp);
set(gca,'YDir','normal')

set(gca,'Xtick',xTick_time+0.5, 'XtickLabel',xTick_time) 
set(gca,'xTickLabel',xTick_label)
set(gca,'Ytick',[1:N_grp])

%% define the color map
color_map = [1 1 1;0.9290, 0.6940, 0.1250];
colormap(color_map)

hold on
%% add color bar edge lines
edge_x = repmat((0:T_otg)+0.5,N_grp+1,1);
edge_y = transpose( repmat((0:N_grp)+0.5,T_otg+1,1) );

plot(edge_x,edge_y,'color',[0.5,0.5,0.5],'LineWidth',1); % vertical lines
hold on
plot(edge_x.',edge_y.','color',[0.5,0.5,0.5],'LineWidth',1); % horizontal lines
hold off

xlabel('Time (hour)')
ylabel('Group Index')

set(gcf,'Color',[1,1,1])
set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',1) set the axis & title & legend font
title(['Served Group',' @ ',mark_title])
saveas(S1_fg_Ug,strcat(output_fd, '\','served_group','_',mark_case,'.fig'));
