% Clear everything to ensure an empty environment before running 
clear all
close all
clc

% Load the Calcium Model's data
load('Sim_TS_Less_Tunnelling.mat')

%% Specify the mesh points to track the calcium concentration in the model
% The tracking points can be changed by referring to the generated mesh in
% the Geom_Mesh code.

% Calcium conc in the cytoplasm
% On top of cortical IP3R
cortical_Cyto_1_p_pointer = Cyto_p(1,:) >= 1443.5 & Cyto_p(1,:) <= 1445 &...
    Cyto_p(2,:) >= 852 & Cyto_p(2,:) <= 853;
cortical_Cyto_1_p = Cyto_p(3,cortical_Cyto_1_p_pointer);
cortical_Cyto_1 = store_res_Cyto_cal(:,cortical_Cyto_1_p);

% On top of deep IP3R 3
deep_Cyto_1_p_pointer = Cyto_p(1,:) >= 1419 & Cyto_p(1,:) <= 1420 &...
    Cyto_p(2,:) >= 150 & Cyto_p(2,:) <= 151;
deep_Cyto_1_p = Cyto_p(3,deep_Cyto_1_p_pointer);
deep_Cyto_1 = store_res_Cyto_cal(:,deep_Cyto_1_p);

%Calcium conc in the ER
% Under the cortical STIM-ORAI and SERCA
cortical_ER_1_p_pointer = ER_p(1,:) >= 444 & ER_p(1,:) <= 445 &...
    ER_p(2,:) >= 874 & ER_p(2,:) <= 875;
cortical_ER_1_p = ER_p(3,cortical_ER_1_p_pointer);
cortical_ER_1 = store_res_ER_cal(:,cortical_ER_1_p);

% Below the cortical IP3R
cortical_ER_2_p_pointer = ER_p(1,:) >= 1452 & ER_p(1,:) <= 1453 &...
    ER_p(2,:) >= 784 & ER_p(2,:) <= 785;
cortical_ER_2_p = ER_p(3,cortical_ER_2_p_pointer);
cortical_ER_2 = store_res_ER_cal(:,cortical_ER_2_p);

% Under deep region IP3R 1
deep_ER_1_p_pointer = ER_p(1,:) >= 1143 & ER_p(1,:) <= 1144 &...
    ER_p(2,:) >= 267 & ER_p(2,:) <= 268;
deep_ER_1_p = ER_p(3,deep_ER_1_p_pointer);
deep_ER_1 = store_res_ER_cal(:,deep_ER_1_p);

% Specify the simulation time step
plot_time = zeros(num_t+1,1);
plot_time(1) = 0;
plot_time(2:num_t+1) = time;

% Write the calcium traces at each point into csv files
cal_point_traces = [plot_time cortical_Cyto_1 deep_Cyto_1 cortical_ER_1 cortical_ER_2 deep_ER_1];
cal_point_traces = cal_point_traces(1:40:end, :);

writematrix(cal_point_traces, 'Fig_3_less_tunnelling_Ca_traces.csv')

%% Plot the time series of the calcium concentration at the specified points
figure(1)
plot(plot_time(1:size(store_res_Cyto_cal,1))/60,...
            cortical_Cyto_1(1:size(store_res_Cyto_cal,1)),'Color','black','LineWidth',3)
hold on
plot(plot_time(1:size(store_res_Cyto_cal,1))/60,...
            deep_Cyto_1(1:size(store_res_Cyto_cal,1)),'Color','red','LineWidth',2)
hold on
line([0,10],[3.1,3.1],'Color','red','LineWidth',3)
hold on
line([0,11],[3.8,3.8],'Color','blue','LineWidth',3)
hold on
line([11,11.5],[3.1,3.1],'Color','black','LineWidth',3)
hold on
line([18,18.5],[3.1,3.1],'Color','black','LineWidth',3)
hold on
line([11,25],[3.8,3.8],'Color',[0.9290 0.6940 0.1250],'LineWidth',3)
text(12,3.1,'His','FontSize',14)
text(19,3.1,'His','FontSize',14)
text(5.5,4.1,'0 Ca^2^+','Color','blue','FontSize',14)
text(5,3.4,'CPA','Color','red','FontSize',14)
text(18,4.1,'Ca^2^+','Color',[0.9290 0.6940 0.1250],'FontSize',14)
xlabel('Time (min)')
ylabel('Ca^2^+ Conc (\muM)')
xlim([0 25])
ylim([0 4.5])
title('Cytosolic Ca^2^+ Traces','FontSize',16)
legend('C1','C2','location','eastoutside')
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
