% Clear everything to ensure an empty environment before running 
clear all
close all
clc

% Load the Calcium Model's data
load('Sim_TS_No_Tunnelling.mat')

%% Specify the mesh points to track the calcium concentration in the model
% The tracking points can be changed by referring to the generated mesh in
% the Geom_Mesh code.

% Calcium conc in the cytoplasm
% Within the cytosolic SOCE microdomain
cortical_Cyto_SOCE_p_pointer = Cyto_p(1,:) >= 600 & Cyto_p(1,:) <= 601 &...
    Cyto_p(2,:) >= 893 & Cyto_p(2,:) <= 894;
cortical_Cyto_SOCE_p = Cyto_p(3,cortical_Cyto_SOCE_p_pointer);
cortical_Cyto_SOCE = store_res_Cyto_cal(:,cortical_Cyto_SOCE_p);

% On top of cortical IP3R
cortical_Cyto_1_p_pointer = Cyto_p(1,:) >= 1502 & Cyto_p(1,:) <= 1503 &...
    Cyto_p(2,:) >= 826 & Cyto_p(2,:) <= 827;
cortical_Cyto_1_p = Cyto_p(3,cortical_Cyto_1_p_pointer);
cortical_Cyto_1 = store_res_Cyto_cal(:,cortical_Cyto_1_p);

% On top of deep IP3R 3
deep_Cyto_1_p_pointer = Cyto_p(1,:) >= 1201 & Cyto_p(1,:) <= 1202 &...
    Cyto_p(2,:) >= 150 & Cyto_p(2,:) <= 151;
deep_Cyto_1_p = Cyto_p(3,deep_Cyto_1_p_pointer);
deep_Cyto_1 = store_res_Cyto_cal(:,deep_Cyto_1_p);

%Calcium conc in the ER
% Under the cortical STIM-ORAI and SERCA
cortical_ER_1_p_pointer = ER_p(1,:) >= 600 & ER_p(1,:) <= 601 &...
    ER_p(2,:) >= 874 & ER_p(2,:) <= 875;
cortical_ER_1_p = ER_p(3,cortical_ER_1_p_pointer);
cortical_ER_1 = store_res_ER_cal(:,cortical_ER_1_p);

% Below the cortical IP3R
cortical_ER_2_p_pointer = ER_p(1,:) >= 1493 & ER_p(1,:) <= 1494 &...
    ER_p(2,:) >= 789 & ER_p(2,:) <= 790;
cortical_ER_2_p = ER_p(3,cortical_ER_2_p_pointer);
cortical_ER_2 = store_res_ER_cal(:,cortical_ER_2_p);

% Under deep region IP3R 1
deep_ER_1_p_pointer = ER_p(1,:) >= 1189 & ER_p(1,:) <= 1190 &...
    ER_p(2,:) >= 274 & ER_p(2,:) <= 275;
deep_ER_1_p = ER_p(3,deep_ER_1_p_pointer);
deep_ER_1 = store_res_ER_cal(:,deep_ER_1_p);

% Specify the simulation time step
plot_time = zeros(num_t+1,1);
plot_time(1) = 0;
plot_time(2:num_t+1) = time;

% Write the calcium traces at each point into csv files
cal_point_traces = [plot_time cortical_Cyto_SOCE cortical_Cyto_1 deep_Cyto_1 cortical_ER_1 cortical_ER_2 deep_ER_1];
cal_point_traces = cal_point_traces(1:40:end, :);

writematrix(cal_point_traces, 'Fig_1_no_tunnelling_Ca_traces.csv')

%% Plot the time series of the calcium concentration at the specified points
figure(1)
plot(plot_time(1:size(store_res_Cyto_cal,1))/60,...
            cortical_Cyto_SOCE(1:size(store_res_Cyto_cal,1)),'Color','black','LineWidth',3)
hold on
plot(plot_time(1:size(store_res_Cyto_cal,1))/60,...
            cortical_Cyto_1(1:size(store_res_Cyto_cal,1)),'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
hold on
plot(plot_time(1:size(store_res_Cyto_cal,1))/60,...
            deep_Cyto_1(1:size(store_res_Cyto_cal,1)),'Color','red','LineWidth',2)
line([0,10],[5.6,5.6],'Color','red','LineWidth',3)
hold on
line([0,11],[6.3,6.3],'Color','blue','LineWidth',3)
hold on
line([11,11.5],[5.6,5.6],'Color','black','LineWidth',3)
hold on
line([18,18.5],[5.6,5.6],'Color','black','LineWidth',3)
hold on
line([11,25],[6.3,6.3],'Color',[0.9290 0.6940 0.1250],'LineWidth',3)
text(12,5.6,'His','FontSize',14)
text(19,5.6,'His','FontSize',14)
text(5.5,6.6,'0 Ca^2^+','Color','blue','FontSize',14)
text(5,5.9,'CPA','Color','red','FontSize',14)
text(18,6.6,'Ca^2^+','Color',[0.9290 0.6940 0.1250],'FontSize',14)
xlabel('Time (min)')
ylabel('Ca^2^+ Conc (\muM)')
xlim([0 25])
ylim([0 7])
title('Cytosolic Ca^2^+ Traces','FontSize',16)
legend('C1','C2','C3','location','eastoutside')
lgd = legend;
lgd.FontSize = 15;

figure(2)
plot(plot_time(1:size(store_res_Cyto_cal,1))/60,...
            cortical_ER_1(1:size(store_res_Cyto_cal,1)),'Color','black','LineWidth',3)
hold on
plot(plot_time(1:size(store_res_Cyto_cal,1))/60,...
            cortical_ER_2(1:size(store_res_Cyto_cal,1)),'Color',[0.9290 0.6940 0.1250],'LineWidth',3)
hold on
plot(plot_time(1:size(store_res_Cyto_cal,1))/60,...
            deep_ER_1(1:size(store_res_Cyto_cal,1)),'Color','red','LineWidth',2)
hold on
line([0,11],[280,280],'Color','blue','LineWidth',3)
hold on
line([0,10],[260,260],'Color','red','LineWidth',3)
hold on
line([11,11.5],[260,260],'Color','black','LineWidth',3)
hold on
line([18,18.5],[260,260],'Color','black','LineWidth',3)
hold on
line([11,25],[280,280],'Color',[0.9290 0.6940 0.1250],'LineWidth',3)
text(12,260,'His','FontSize',14)
text(19,260,'His','FontSize',14)
text(5.5,290,'0 Ca^2^+','Color','blue','FontSize',14)
text(5,268,'CPA','Color','red','FontSize',14)
text(20,290,'Ca^2^+','Color',[0.9290 0.6940 0.1250],'FontSize',14)
xlabel('Time (min)')
ylabel('Ca^2^+ Conc (\muM)')
xlim([0 25])
ylim([0 300])
title('ER Ca^2^+ Traces','FontSize',16)
legend('E1','E2','E3','location','eastoutside')
lgd = legend;
lgd.FontSize = 15;