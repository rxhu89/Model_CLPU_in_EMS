%% Coded by Rongxing Hu in May 2024, guliusha@outlook.com, NCSU
%% Coded and tested on Matlab 2019
%% compare and plot the deviations of HVAC simulations and EMS results

close all

N_hour = N_step;
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
    fig_size_s = [0.3,0.4,0.3,0.45]; % 
    fig_size_b = [0.3,0.4,0.3,0.5]; % 
    xTick_time = [0:4:T_otg]./1;
    xTick_label = vec_time_label(1+[0:4:T_otg]);
elseif hrzn_d >1 && hrzn_d <=2
    fig_size   = [0.3,0.5,0.4,0.4]; % for 2-3 day data; served node is different
    fig_size_s = [0.3,0.4,0.5,0.45]; % 
    fig_size_b = [0.3,0.4,0.4,0.5]; % 
    xTick_time = [0:4:T_otg]./1;
    xTick_label = vec_time_label(1+[0:4:T_otg]);
elseif hrzn_d >=3 && hrzn_d <=4
    fig_size   = [0.2,0.5,0.6,0.4]; % for 4-5 day data; served node is different
    fig_size_s = [0.2,0.4,0.7,0.45]; % 
    fig_size_b = [0.2,0.4,0.6,0.5]; % 
    xTick_time = [0:6:T_otg]./1;
    xTick_label = vec_time_label(1+[0:6:T_otg]);
else
    fig_size   = [0.2,0.5,0.7,0.4]; % for 6-7 day data; served node is different
    fig_size_s = [0.1,0.4,0.8,0.45]; % 
    fig_size_b = [0.1,0.4,0.7,0.5]; % 
    xTick_time = [0:6:T_otg]./1;
    xTick_label = vec_time_label(1+[0:6:T_otg]);
end

%% plot power sum
simu_Phvac_1 = sum(house_pow,2)';
simu_Phvac_60 = sum(hvac_house_pow_60_ori,1);
EMS_Phvac_60 = sum(sol_Pserve_hvac,1);

% house_grp_ph
map_hs_grp_ph = zeros(N_house,3); % [hs_idx, grp_idx, ph_idx]
temp_order_start = 1;
for i = 1:N_grp
    for j = 1:3
        temp_house_idx_ph_grp = M_house_idx_my_grp_ph((i-1)*3+j,:);
        temp_house_idx_ph_grp(temp_house_idx_ph_grp==0) = []; % delete zeros
        temp_N_hs_in_grp_ph = length(temp_house_idx_ph_grp);
        
        map_hs_grp_ph(temp_order_start:temp_order_start+temp_N_hs_in_grp_ph-1,:) = [temp_house_idx_ph_grp', i*ones(temp_N_hs_in_grp_ph,1), j*ones(temp_N_hs_in_grp_ph,1)] ;
        temp_order_start = temp_order_start + temp_N_hs_in_grp_ph;
    end
end
map_hs_grp_ph = sortrows(map_hs_grp_ph,1);

fg_pow_all = figure;
set(gcf,'unit','normalized','position',fig_size); % 
yyaxis left
bar([1:1:N_hour]-0.5,EMS_Phvac_60,1,'facecolor',[0.3010 0.7450 0.9330],'edgecolor',[0.3010 0.7450 0.9330]) %,'-m','LineWidth',1)
hold on
plot([1:1:N_hour*60]./60,simu_Phvac_1,'-r','LineWidth',1)
hold on
plot([1:1:N_hour],simu_Phvac_60,'LineWidth',1.2) %,'-m','LineWidth',1)

ylim([0,3500])
ylabel('HVAC Power(kW)')

yyaxis right
plot([1:N_hour]-0.5,Tout_data_hour,'-->','LineWidth',1.5,'color',[1,0.7,0.0] ) %orange [0.85,0.33,0.10]
ylabel(strcat('Temperature  (', char(0176),'C',')') )
ylim([20,40])

title(['Simulated HVACs',' @ ', mark_case_title_simu ])
legend('Schedule  (60-min)','Simu (1-min)','Simu (60-min)','Temperature','Location','northwest','Orientation','vertical') %'horizontal'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'xTick',xTick_time) % [0:4:T_otg]./2)
set(gca,'xTickLabel',xTick_label)
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14); %,'FontWeight','Bold', 'LineWidth', 2);

saveas(fg_pow_all,strcat(output_fd_hvac_simu, 'HVAC_pow_all','_',mark_case_simu,'.fig'));

%% ===============================================

%% set the components
P_PV_rate = EMS_result.sys_component.P_PV_rate;
P_bat_rate = EMS_result.sys_component.P_bat_rate;
E_bat_rate = EMS_result.sys_component.E_bat_rate;
k_dc   = EMS_result.sys_component.k_dc;
soc_0  = EMS_result.opti_ini.soc_0;

soc_ini = soc_0;

input_fd_hvac_norm = strcat(input_fd_hvac,'Output_HVAC_normal','\');

%% Process data
simu_Phvac_1 = sum(house_pow,2)';          % hvac_house_pow_60_ori [N_step, N_house]
simu_Phvac_60 = sum(hvac_house_pow_60_ori,1); % hvac_house_pow_60_ori [N_house,N_step]
EMS_Phvac_60 = sum(sol_Pserve_hvac,1);     % sol_Pserve_hvac [N_grp*3,N_step]

sol_Pclpu_grp = EMS_result.sol_Pclpu_grp;    % this is the group-ph load

% ==== select data length, day, hour
day1  = d ;
hour1 = otg_sh ;

% -- import normal HVAC load
norm_Phvac_hs_1 = importdata( strcat(input_fd_hvac_norm,'bus123_hs_Phvac_Aug_m_1_rand.mat') );
norm_Phvac_hs_1 = norm_Phvac_hs_1( ((day1-1)*24+hour1-1)*60+1:((day1-1)*24+hour1-1+hrzn_h)*60, : ); % [N_step, N_hs]
norm_Phvac_hs_1 = norm_Phvac_hs_1 .* sol_Uhouse_minu';

% -- to group data
simu_Phvac_grp_ph_1 = zeros(3*N_grp, size(house_pow,1));
simu_Phvac_grp_ph_60 = zeros(3*N_grp, size(hvac_house_pow_60_ori,2));
norm_Phvac_grp_ph_1 = zeros(3*N_grp, size(norm_Phvac_hs_1,1));
for i = 1:N_grp
    for j = 1:3
        temp_house_idx_ph_grp = M_house_idx_my_grp_ph((i-1)*3+j,:);
        temp_house_idx_ph_grp(temp_house_idx_ph_grp==0) = []; % delete zeros
        
        simu_Phvac_grp_ph_1((i-1)*3+j,:) = sum(house_pow(:,temp_house_idx_ph_grp),2)' ;
        simu_Phvac_grp_ph_60((i-1)*3+j,:) = sum(hvac_house_pow_60_ori(temp_house_idx_ph_grp,:));
        norm_Phvac_grp_ph_1((i-1)*3+j,:) = sum(norm_Phvac_hs_1(:,temp_house_idx_ph_grp),2)' ;
    end
end
simu_Pclpu_grp_ph_1 = simu_Phvac_grp_ph_1 - norm_Phvac_grp_ph_1;

simu_Phvac_pcc_1 = [sum(simu_Phvac_grp_ph_1(1:3:end,:)); sum(simu_Phvac_grp_ph_1(2:3:end,:)); sum(simu_Phvac_grp_ph_1(3:3:end,:)) ];
simu_Phvac_pcc_60 = [sum(simu_Phvac_grp_ph_60(1:3:end,:)); sum(simu_Phvac_grp_ph_60(2:3:end,:)); sum(simu_Phvac_grp_ph_60(3:3:end,:)) ];
EMS_Phvac_pcc_60 = [sum(sol_Pserve_hvac(1:3:end,:)); sum(sol_Pserve_hvac(2:3:end,:)); sum(sol_Pserve_hvac(3:3:end,:)) ];
norm_Phvac_pcc_1 = [sum(norm_Phvac_grp_ph_1(1:3:end,:)); sum(norm_Phvac_grp_ph_1(2:3:end,:)); sum(norm_Phvac_grp_ph_1(3:3:end,:)) ];

simu_Phvac_grp_1  = zeros(N_grp, size(simu_Phvac_pcc_1,2));
simu_Phvac_grp_60 = zeros(N_grp, size(simu_Phvac_pcc_60,2));

EMS_Phvac_grp_60  = zeros(N_grp, size(EMS_Phvac_pcc_60,2));
EMS_Pclpu_grp_60  = zeros(N_grp, size(EMS_Phvac_pcc_60,2));

norm_Phvac_grp_1  = zeros(N_grp, size(norm_Phvac_pcc_1,2));
simu_Pclpu_grp_1  = zeros(N_grp, size(norm_Phvac_pcc_1,2));
for i = 1:N_grp
    simu_Phvac_grp_1(i,:) = sum( simu_Phvac_grp_ph_1((i-1)*3+1:i*3,:) );
    simu_Phvac_grp_60(i,:) = sum( simu_Phvac_grp_ph_60((i-1)*3+1:i*3,:) );
    
    EMS_Phvac_grp_60(i,:) = sum( sol_Pserve_hvac((i-1)*3+1:i*3,:) );
    EMS_Pclpu_grp_60(i,:) = sum( sol_Pclpu_grp((i-1)*3+1:i*3,:) );
    
    norm_Phvac_grp_1(i,:) = sum( norm_Phvac_grp_ph_1((i-1)*3+1:i*3,:) );
    simu_Pclpu_grp_1(i,:) = sum( simu_Pclpu_grp_ph_1((i-1)*3+1:i*3,:) );
end

%% plot Phvac of groups for comparison
fg_Phavc_grp = figure;
set(gcf,'unit','normalized','position',fig_size_b); % the size of the figure

sgtitle(['Power of HVACs',' @ ', mark_case_title_simu ])
subplot(2,1,1)
yyaxis left
hB = bar([1:1:N_hour]-0.5,EMS_Phvac_grp_60,1,'stacked','FaceAlpha',1); %,'-m','LineWidth',1) 'edgecolor','None','facecolor',[0.3010 0.7450 0.9330],
color_order = lines(N_grp);
for g = 1:N_grp
    hB(g).FaceColor = color_order(g,:);
end
ylabel(strcat('HVAC Power (kW)') )
ylim([0,3500])
xlim([0,N_hour])

yyaxis right
plot([1:N_hour]-0.5,Tout_data_hour,'-->','LineWidth',1.5,'color',[1,0.7,0.0] ) %orange [0.85,0.33,0.10]
ylabel(strcat('Temperature  (', char(0176),'C',')') )
ylim([20,40])

legend('G1 (EMS)','G2 (EMS)','G3 (EMS)','G4 (EMS)','G5 (EMS)','Temperature','Location','northwest','Orientation','horizontal','NumColumns', 5) %'horizontal' 'vertical'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'xTick',xTick_time) % [0:4:T_otg]./2)
set(gca,'xTickLabel',xTick_label)
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14); %,'FontWeight','Bold', 'LineWidth', 2);

subplot(2,1,2)
plot([1:1:N_hour*60]./60,simu_Phvac_grp_1,'-','LineWidth',1) %,'-m','LineWidth',1)
ylabel(strcat('HVAC Power (kW)') )
ylim([-200,1500])
xlim([0,N_hour])

legend('G1 (Simu)','G2 (Simu)','G3 (Simu)','G4 (Simu)','G5 (Simu)','Location','northwest','Orientation','horizontal','NumColumns', 5) %'horizontal' 'vertical'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'xTick',xTick_time) % [0:4:T_otg]./2)
set(gca,'xTickLabel',xTick_label)
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14); %,'FontWeight','Bold', 'LineWidth', 2);

saveas(fg_Phavc_grp,strcat(output_fd_hvac_simu, 'Phvac_grp','_',mark_case_simu,'.fig'));

%% plot Pclpu of groups for comparison
fg_Phavc_clpu_grp = figure;
set(gcf,'unit','normalized','position',fig_size_b); % the size of the figure

sgtitle(['Power of CLPU',' @ ', mark_case_title_simu ])
subplot(2,1,1)
yyaxis left
hB = bar([1:1:N_hour]-0.5,EMS_Pclpu_grp_60,1,'stacked','FaceAlpha',1); %,'-m','LineWidth',1) 'edgecolor','None','facecolor',[0.3010 0.7450 0.9330],
color_order = lines(N_grp);
for g = 1:N_grp
    hB(g).FaceColor = color_order(g,:);
end
ylabel(strcat('CLPU Power (kW)') )
ylim([0,3500])
xlim([0,N_hour])

yyaxis right
plot([1:N_hour]-0.5,Tout_data_hour,'-->','LineWidth',1.5,'color',[1,0.7,0.0] ) %orange [0.85,0.33,0.10]
ylabel(strcat('Temperature  (', char(0176),'C',')') )
ylim([20,40])

legend('G1 (EMS)','G2 (EMS)','G3 (EMS)','G4 (EMS)','G5 (EMS)','Temperature','Location','northwest','Orientation','horizontal','NumColumns', 5) %'horizontal' 'vertical'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'xTick',xTick_time) % [0:4:T_otg]./2)
set(gca,'xTickLabel',xTick_label)
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14); %,'FontWeight','Bold', 'LineWidth', 2);

subplot(2,1,2)
plot([1:1:N_hour*60]./60,simu_Pclpu_grp_1,'-','LineWidth',1) %,'-m','LineWidth',1)
ylabel(strcat('CLPU Power (kW)') )
ylim([-200,1500])
xlim([0,N_hour])

legend('G1 (Simu)','G2 (Simu)','G3 (Simu)','G4 (Simu)','G5 (Simu)','Location','northwest','Orientation','horizontal','NumColumns', 5) %'horizontal' 'vertical'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'xTick',xTick_time) % [0:4:T_otg]./2)
set(gca,'xTickLabel',xTick_label)
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14); %,'FontWeight','Bold', 'LineWidth', 2);

saveas(fg_Phavc_clpu_grp,strcat(output_fd_hvac_simu, 'Phvac_clpu_grp','_',mark_case_simu,'.fig'));

%% import EMS-hvac results @ BESS
sol_Ebat = EMS_result.sol_Ebat;    % Total
sol_Pc = EMS_result.sol_Pc;    % 3-ph
sol_Pd = EMS_result.sol_Pd;    % 3-ph

%%  Obtain the BESS charging and discharging power baesed on 
EMS_Pbat_60  = sol_Pd - sol_Pc;
EMS_Pbat_1 = kron(EMS_Pbat_60,ones(1,60));
EMS_Phvac_pcc_1 = kron(EMS_Phvac_pcc_60,ones(1,60));

delta_Phvac_pcc_1 = simu_Phvac_pcc_1 - EMS_Phvac_pcc_1;
simu_Pbat_1 = EMS_Pbat_1 + delta_Phvac_pcc_1;

simu_Pbat_ch_1 = simu_Pbat_1;
simu_Pbat_ch_1(simu_Pbat_ch_1>0) = 0;
simu_Pbat_ch_1 = -1 * simu_Pbat_ch_1;

simu_Pbat_dch_1 = simu_Pbat_1;
simu_Pbat_dch_1(simu_Pbat_dch_1 <0) = 0;

delta_Ebat_pcc_1 = ( (1/k_dc)*simu_Pbat_dch_1 - k_dc*simu_Pbat_ch_1 ) / 60; % 3-ph @ To hour @ kWh
delta_Ebat_1 = sum(delta_Ebat_pcc_1,1);
simu_Ebat_1 = soc_ini*E_bat_rate - cumsum(delta_Ebat_1);

%% Save data
Pbat_havc_simu_based.simu_Pbat_1  = simu_Pbat_1;
Pbat_havc_simu_based.simu_Pbat_ch_1  = simu_Pbat_ch_1;
Pbat_havc_simu_based.simu_Pbat_dch_1  = simu_Pbat_dch_1;
Pbat_havc_simu_based.delta_Ebat_pcc_1  = delta_Ebat_pcc_1;
Pbat_havc_simu_based.simu_Ebat_1  = simu_Ebat_1;

save(strcat(output_fd_hvac_simu,  'Pbat_havc_simu_based_1m','_',mark_case_simu,'.mat'), 'Pbat_havc_simu_based');

Phvac_clpu_simu.simu_Phvac_grp_ph_1 = simu_Phvac_grp_ph_1;
Phvac_clpu_simu.simu_Phvac_grp_ph_60 = simu_Phvac_grp_ph_60;
Phvac_clpu_simu.norm_Phvac_grp_ph_1 = norm_Phvac_grp_ph_1;
Phvac_clpu_simu.simu_Pclpu_grp_ph_1 = simu_Pclpu_grp_ph_1;

save(strcat(output_fd_hvac_simu,  'Phvac_clpu_simu','_',mark_case_simu,'.mat'), 'Phvac_clpu_simu');

%% plot bess overall
fg_Phavc_grp = figure;
set(gcf,'unit','normalized','position',fig_size); % the size of the figure

yyaxis left
bar([1:1:N_hour]-0.5,sol_Ebat,1,'facecolor',[0.3010 0.7450 0.9330],'edgecolor','None','FaceAlpha',0.5) %,'-m','LineWidth',1)
hold on
plot([1:1:N_hour*60]./60,simu_Ebat_1,'k-','LineWidth',1) %,'-m','LineWidth',1)
ylabel(strcat('BESS Energy (kWh)') )
ylim([0,E_bat_rate*1.1])

yyaxis right
plot([1:1:N_hour*60]./60,sum(EMS_Pbat_1,1),'b-','LineWidth',1)
hold on
plot([1:1:N_hour*60]./60,sum(simu_Pbat_1,1),':','LineWidth',1.5)
ylim([-P_bat_rate*3-500,P_bat_rate*3+500])
ylabel('BESS Power(kW) (dch+)')

title(['BESS',' @ ', mark_case_title_simu ])
legend('Energy-EMS (60-min)','Energy-Simu (1-min)','pow-EMS (60-min)','pow-Simu (1-min)','Location','southwest','Orientation','vertical') %'horizontal'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'xTick',xTick_time) % [0:4:T_otg]./2)
set(gca,'xTickLabel',xTick_label)
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14); %,'FontWeight','Bold', 'LineWidth', 2);

saveas(fg_Phavc_grp,strcat(output_fd_hvac_simu, 'BESS_overall','_',mark_case_simu,'.fig'));

%% plot bess pow
fg_bat_pow = figure;
set(gcf,'unit','normalized','position',fig_size); % the size of the figure

plot([1:1:N_hour*60]./60,EMS_Pbat_1(1:3:end,:),'-','LineWidth',1,'color',[1.0,0.7,0.0] )
hold on
plot([1:1:N_hour*60]./60,EMS_Pbat_1(2:3:end,:),'g-','LineWidth',1 )
hold on
plot([1:1:N_hour*60]./60,EMS_Pbat_1(3:3:end,:),'r-','LineWidth',1 )
hold on

plot([1:1:N_hour*60]./60,simu_Pbat_1(1:3:end,:),':','LineWidth',1.5,'color',[1.0,0.7,0.0] )
hold on
plot([1:1:N_hour*60]./60,simu_Pbat_1(2:3:end,:),'g:','LineWidth',1.5 )
hold on
plot([1:1:N_hour*60]./60,simu_Pbat_1(3:3:end,:),'r:','LineWidth',1.5 )

ylim([-1*P_bat_rate-500,P_bat_rate+500])
ylabel('BESS Power(kW) (dch+)')

title(['BESS',' @ ', mark_case_title_simu ])
legend('pow-EMS a','pow-EMS b','pow-EMS c','pow-Simu a','pow-Simu b','pow-Simu c','Location','northwest','Orientation','vertical') %'horizontal'
hLegend = legend;
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
set(gca,'xTick',xTick_time) % [0:4:T_otg]./2)
set(gca,'xTickLabel',xTick_label)
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14); %,'FontWeight','Bold', 'LineWidth', 2);

saveas(fg_bat_pow,strcat(output_fd_hvac_simu, 'BESS_pow','_',mark_case_simu,'.fig'));
