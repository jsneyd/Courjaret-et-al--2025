% Clear everything to ensure an empty environment before running 
clear all
close all
clc

% Load the Calcium Model's data
load('Sim_TS_Tunnelling.mat')

% Specify the simulation time step used
plot_time = zeros(num_t+1,1);
plot_time(1) = 0;
plot_time(2:num_t+1) = time;

% Write the video destination and file name for the cytosolic video
v_cyto_pic = VideoWriter("/Users/llee544/Desktop/Tunnel_Cyto","MPEG-4");
v_cyto_pic.FrameRate = 40; % Specify the desired frame rate per second
open(v_cyto_pic); % Open the video file to store the generated video

%Parameter that reduces the video size and length by reducing the amount of
%data used per video (reducing this parameter and increasing the FrameRate
%produces a high definition video but at the cost of memory and computation
%time)
k_1 = 200; 

% Generate a video that tracks the calcium conc in the cytosolic domain
fig_1 = figure('position',[50,50,960,900]);
for i = 1:(size(store_res_Cyto_cal,1)-1)/k_1+1
    if i == 1
        pdeplot(Cyto_model.Mesh,XYData = store_res_Cyto_cal(1,:),ColorMap='jet')
        xlim([-50 2400])
        ylim([-50 950])
        clim([0 4])
        set(gca,'YTick', [],'XTick',[])
    elseif i == (size(store_res_Cyto_cal,1)-1)/k_1+1
        pdeplot(Cyto_model.Mesh,XYData = store_res_Cyto_cal(size(store_res_Cyto_cal,1),:),ColorMap='jet')
        xlim([-50 2400])
        ylim([-50 950])
        clim([0 4])
        set(gca,'YTick', [],'XTick',[])
    else
        pdeplot(Cyto_model.Mesh,XYData = store_res_Cyto_cal(i*k_1,:),ColorMap='jet')
        xlim([-50 2400])
        ylim([-50 950])
        clim([0 4])
        set(gca,'YTick', [],'XTick',[])
    end
    
    if i == 1
        sgtitle("Cytoplasm, Time = 0 min (CPA On, Ext Ca^2^+ Off, His Off)","FontSize",20)
    elseif i <= 600/t_step/k_1
        sgtitle("Cytoplasm, Time = "+round((i-1)*k_1*t_step/60,1)+" min (CPA On, Ext Ca^2^+ Off, His Off)","FontSize",20)
    elseif i > 600/t_step/k_1 && i < 660/t_step/k_1
        sgtitle("Cytoplasm, Time = "+round((i-1)*k_1*t_step/60,1)+" min (CPA Off, Ext Ca^2^+ Off, His Off)","FontSize",20)
    elseif (i >= 660/t_step/k_1 && i <= 690/t_step/k_1) ||...
           (i >= 1080/t_step/k_1 && i <= 1110/t_step/k_1)
        sgtitle("Cytoplasm, Time = "+round((i-1)*k_1*t_step/60,1)+" min (CPA Off, Ext Ca^2^+ On, His On)","FontSize",20)
    elseif i == (size(store_res_Cyto_cal,1)-1)/k_1+1
        sgtitle("Cytoplasm, Time = "+round((i-1)*k_1*t_step/60,1)+" min (CPA Off, Ext Ca^2^+ On, His Off)","FontSize",20)
    else
        sgtitle("Cytoplasm, Time = "+round((i-1)*k_1*t_step/60,1)+" min (CPA Off, Ext Ca^2^+ On, His Off)","FontSize",20)
    end
    
    f_1(i) = getframe(fig_1);
    clf
end

writeVideo(v_cyto_pic,f_1) % Store the video in the specified file and destination

close(v_cyto_pic); % Close the video file after generating the video

% Write the video destination and file name for the ER video
v_ER_pic = VideoWriter("/Users/llee544/Desktop/Tunnel_ER","MPEG-4");
v_ER_pic.FrameRate = 40; % Specify the desired frame rate per second
open(v_ER_pic); % Open the video file to store the generated video

% Generate a video that tracks the calcium conc in the ER domain
fig_2 = figure('position',[50,50,960,900]);
for i = 1:(size(store_res_Cyto_cal,1)-1)/k_1+1
    if i == 1
        pdeplot(ER_model.Mesh,XYData = store_res_ER_cal(1,:),ColorMap='jet')
        xlim([-50 2400])
        ylim([-50 950])
        clim([0 220])
        set(gca,'YTick', [],'XTick',[])
    elseif i == (size(store_res_Cyto_cal,1)-1)/k_1+1
        pdeplot(ER_model.Mesh,XYData = store_res_ER_cal(size(store_res_Cyto_cal,1),:),ColorMap='jet')
        xlim([-50 2400])
        ylim([-50 950])
        clim([0 220])
        set(gca,'YTick', [],'XTick',[])
    else
        pdeplot(ER_model.Mesh,XYData = store_res_ER_cal(i*k_1,:),ColorMap='jet')
        xlim([-50 2400])
        ylim([-50 950])
        clim([0 220])
        set(gca,'YTick', [],'XTick',[])
    end
    
    if i == 1
        sgtitle("ER, Time = 0 min (CPA On, Ext Ca^2^+ Off, His Off)","FontSize",20)
    elseif i <= 600/t_step/k_1
        sgtitle("ER, Time = "+round((i-1)*k_1*t_step/60,1)+" min (CPA On, Ext Ca^2^+ Off, His Off)","FontSize",20)
    elseif i > 600/t_step/k_1 && i < 660/t_step/k_1
        sgtitle("ER, Time = "+round((i-1)*k_1*t_step/60,1)+" min (CPA Off, Ext Ca^2^+ Off, His Off)","FontSize",20)
    elseif (i >= 660/t_step/k_1 && i <= 690/t_step/k_1) ||...
           (i >= 1080/t_step/k_1 && i <= 1110/t_step/k_1)
        sgtitle("ER, Time = "+round((i-1)*k_1*t_step/60,1)+" min (CPA Off, Ext Ca^2^+ On, His On)","FontSize",20)
    elseif i == (size(store_res_Cyto_cal,1)-1)/k_1+1
        sgtitle("ER, Time = "+round((i-1)*k_1*t_step/60,1)+" min (CPA Off, Ext Ca^2^+ On, His Off)","FontSize",20)
    else
        sgtitle("ER, Time = "+round((i-1)*k_1*t_step/60,1)+" min (CPA Off, Ext Ca^2^+ On, His Off)","FontSize",20)
    end
    
    f_2(i) = getframe(fig_2);
    clf
end

writeVideo(v_ER_pic,f_2) % Store the video in the specified file and destination

close(v_ER_pic); % Close the video file after generating the video

