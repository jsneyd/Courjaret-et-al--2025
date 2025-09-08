%% FEM Code for the No Tunnelling Calcium Model

% Clear everything to ensure an empty environment before running 
clear all
close all
clc

% Track the run time for the simulation
tic

% Load the No Tunnelling Calcium Model's mesh data
load("Mesh_no_tunnel.mat")
% [ER_p,ER_e,ER_t], [Cyto_p,Cyto_e,Cyto_t], inter_bound_p

%% Specify the simulation time parameters
t_start = 0; %Start time
t_end = 1500; %End time
t_step = 0.0025; %Simulation time step size (0.0025 s)
num_t = round(t_end/t_step,0); %Number of time step
time = linspace(t_start+t_step,t_end,num_t); %Vector that contains each simulation time point

%% Initial Conditions
%Initial conditions for h variable in the cortical region
cortical_h_pointer = (Cyto_p(1,:) >= 1300 & Cyto_p(1,:) <= 1700 & Cyto_p(2,:) == 800);
cortical_h_p = round(Cyto_p(:,cortical_h_pointer),8);
cortical_h = zeros(num_t+1,size(cortical_h_p,2));
cortical_h(1,:) = 0.8;
cortical_h_init = cortical_h(1,:)';

%Initial conditions for h variable in the deep region
deep_h_pointer = ((Cyto_p(1,:) >= 550 & Cyto_p(1,:) <= 1850 & Cyto_p(2,:) == 200) |...
    (Cyto_p(1,:) >= 550 & Cyto_p(1,:) <= 1850 & Cyto_p(2,:) == 350) |...
    (Cyto_p(1,:) == 1850 & Cyto_p(2,:) >= 200 & Cyto_p(2,:) <= 350));
deep_h_p = round(Cyto_p(:,deep_h_pointer),8);
deep_h = zeros(num_t+1,size(deep_h_p,2));
deep_h(1,:) = 0.8;
deep_h_init = deep_h(1,:)';

%Initial conditions for c and ce variables in the Cytoplasm and ER region
init_val_ER_cal = init_condition(ER_p,"ER_cal");
init_val_Cyto_cal = init_condition(Cyto_p,"Cyto_cal");

store_res_ER_cal = zeros(num_t+1,size(ER_p,2));
store_res_ER_cal(1,:) = init_val_ER_cal';
store_res_Cyto_cal = zeros(num_t+1,size(Cyto_p,2));
store_res_Cyto_cal(1,:) = init_val_Cyto_cal';

M_ER = MassAssembler2D(ER_p,ER_t);
M_Cyto = MassAssembler2D(Cyto_p,Cyto_t);

A_ER = StiffnessAssembler2D(ER_p,ER_t,@Diff_val,"ER_cal");
A_Cyto = StiffnessAssembler2D(Cyto_p,Cyto_t,@Diff_val,"Cyto_cal");

R_ER = RobinStiff2D(ER_p,ER_e,@kappa_val,"ER_cal");
R_Cyto = RobinStiff2D(Cyto_p,Cyto_e,@kappa_val,"Cyto_cal");

L_ER = LoadAssembler2D(ER_p,ER_t,@f_val,"ER_cal"); % Get the Load vector in subdomain 1
L_Cyto = LoadAssembler2D(Cyto_p,Cyto_t,@f_val,"Cyto_cal"); % Get the Load vector in subdomain 2

%% Solving the PDE system at every time step
% Loop over time
ER_cal_new = init_val_ER_cal;
Cyto_cal_new = init_val_Cyto_cal;
cortical_h_new = cortical_h_init;
deep_h_new = deep_h_init;

for i = 1:num_t
    % Set t-1 solution into a temporary vector
    ER_cal_old = ER_cal_new;
    Cyto_cal_old = Cyto_cal_new;
    cortical_h_old = cortical_h_new;
    deep_h_old = deep_h_new;

    r_ER = RobinLoad2D(ER_p,ER_e,@kappa_val,@gD_val,@gN_val,...
        ER_cal_old,Cyto_cal_old,inter_bound_p,"ER_cal",...
        cortical_h_p,cortical_h_old,deep_h_p,deep_h_old,i,t_step);
    
    r_Cyto = RobinLoad2D(Cyto_p,Cyto_e,@kappa_val,@gD_val,@gN_val,...
        ER_cal_old,Cyto_cal_old,inter_bound_p,"Cyto_cal",...
        cortical_h_p,cortical_h_old,deep_h_p,deep_h_old,i,t_step);
    
    ER_cal_new = (M_ER + t_step*(A_ER+R_ER))\(M_ER*ER_cal_old + t_step*(L_ER+r_ER));
    Cyto_cal_new = (M_Cyto + t_step*(A_Cyto+R_Cyto))\(M_Cyto*Cyto_cal_old + t_step*(L_Cyto+r_Cyto));
    
    cortical_h_new = calculate_h_new(cortical_h_p,Cyto_cal_old,cortical_h_old,t_step);
    deep_h_new = calculate_h_new(deep_h_p,Cyto_cal_old,deep_h_old,t_step);
    
    %Store the result
    store_res_ER_cal(i+1,:) = round(ER_cal_new',8);
    store_res_Cyto_cal(i+1,:) = round(Cyto_cal_new',8);
    
    cortical_h(i+1,:) = cortical_h_new';
    deep_h(i+1,:) = deep_h_new';
    
    i % Print the loop number to ensure the simulation is running.
end

save("Sim_TS_No_Tunnelling.mat",'-v7.3')

check_neg_Cyto = store_res_Cyto_cal < 0;
check_neg_ER = store_res_ER_cal < 0;

% Print any negative output values if there is any 
% (due to simulation time step being too large). Choose a smaller time
% step for simulation if there is negative value.
store_res_Cyto_cal(check_neg_Cyto)
store_res_ER_cal(check_neg_ER)

toc

%% Functions which can be changed to reflect the PDE system
%The initial condition function for c and ce
function z = init_condition(p,compartment)
    nz = size(p,2); %A handy variable that store the number of nodal points in the Point matrix
    z = zeros(nz,1); %Initialize the variable that stores the output
    for i = 1:nz
        if compartment == "ER_cal"
            z(i) = 200; 
        elseif compartment == "Cyto_cal"
            z(i) = 0.1; 
        else
            fprintf('Something is wrong with the Initial Condition function');
            break
        end
    end
end

% The Diffusion coefficient function for the system
function z = Diff_val(x,y,compartment)
    b_Cyto_total = 2000;
    K_Cyto = 9;
    b_ER_total = 2000;
    K_ER = 36;
    D_cal = 223e6;
    
    if compartment == "Cyto_cal"
        z = D_cal*K_Cyto/(K_Cyto+b_Cyto_total);
    elseif compartment == "ER_cal"
        z = D_cal*K_ER/(K_ER+b_ER_total);
    else
        fprintf('Something is wrong with the Diffusion function');
    end
end

% RHS of the PDE system which does not involves the spatial and time derivatives 
function z = f_val(x,y,compartment)
    z = 0;
end

% The Robin parameter to determine whether it is Dirichlet or Neumann
% boundary in the cytoplasm for the specified locations.
function z = kappa_val(x,y,compartment)
    z = 0;
end

% Specify the Dirichlet boundary condition for the cytoplasm domain
function z = gD_val(x,y,compartment)
    z = 0;
end

% Function to get the local calcium at the shared interior boundaries
function [cal_ER,cal_Cyto] = loc_cal_func(inter_bound_p,loc2glb,ER_cal,Cyto_cal,compartment)
    if compartment == "ER_cal"
        compart_point = 4;
    elseif compartment == "Cyto_cal"
        compart_point = 3;
    else
        fprintf('Something is wrong with the local calcium function');
    end
    
    glb_pointer_1 = (inter_bound_p(compart_point,:) == loc2glb(1,:));
    glb_pointer_2 = (inter_bound_p(compart_point,:) == loc2glb(2,:));
    glb_point_nodal_1 = inter_bound_p(:,glb_pointer_1);
    glb_point_nodal_2 = inter_bound_p(:,glb_pointer_2);
    loc_ER_cal_1 = ER_cal(glb_point_nodal_1(4,:));
    loc_ER_cal_2 = ER_cal(glb_point_nodal_2(4,:));
    loc_Cyto_cal_1 = Cyto_cal(glb_point_nodal_1(3,:));
    loc_Cyto_cal_2 = Cyto_cal(glb_point_nodal_2(3,:));

    cal_ER = (loc_ER_cal_1+loc_ER_cal_2)/2;
    cal_Cyto = (loc_Cyto_cal_1+loc_Cyto_cal_2)/2;
end

% Specify the Neumann boundary condition for the ER domain (also
% tracks fluxes in the domain)
function z = gN_val(x,y,loc2glb,...
            ER_cal,Cyto_cal,inter_bound_p,compartment, ...
            cortical_h_p,cortical_h,deep_h_p,deep_h,loop_num,t_step)
    b_Cyto_total = 2000;
    K_Cyto = 9;
    b_ER_total = 2000;
    K_ER = 36;
    
    if compartment == "Cyto_cal"
        buff_scale = K_Cyto/(K_Cyto+b_Cyto_total);
    elseif compartment == "ER_cal"
        buff_scale = K_ER/(K_ER+b_ER_total);
    end
    gamma = (1830500/4056)/(284500/1031);
    %Jserca parameters
    Vp = 240*1667.67;
    K_bar_serca = 2.5e-7;
    k_serca = 0.2;
    %Jpm parameters
    Vpm = 12*4500;
    Kpm = 0.2;
    %Jin parameters
    K_stim = 84.5; 
    alpha0 = 11*850; 
    alpha1 = 12*97280;
    STIM_Hill = 4.2; 
    %Jipr parameters
    kf = 48*442.584; 
    Kc = 0.2;
    Kh = 0.8; 
    Kp = 0.2;
    %SERCA Mapper blocker
    blocker = 1; % Cortical SERCA pump blocker. Set between 0 and 1 to indicate the percentage blocked.
    %Thapsigargin
    thap = 1; % 0 when thap is on and 1 when thap is off
    %Agonist stimulation time period
    loop_time = loop_num*t_step;
    if (loop_time >= 660 && loop_time <= 690) || ...
            (loop_time >= 1080 && loop_time <= 1110)
        Ago = 100; 
    else
        Ago = 0;
    end
    %Allow external calcium period
    if (loop_time >= 660 && loop_time <= 10300)
        ext_cal_on = 1;
    else
        ext_cal_on = 0;
        
    end
    %Block SERCA pump period
    if loop_time >= 0 && loop_time <= 600
        CPA_on = 0; %This means SERCA is blocked
    else
        CPA_on = 1;
    end
    
    % Calculate calcium flux at each IP3R, SERCA, SOCE, and PMCA location
    % Cortical SERCA flux
    if ((x(1) == 550 && y(1) >= 700 && y(1) < 885) && ...
        (x(2) == 550 && y(2) >= 700 && y(2) < 885)) || ...
      ((x(1) == 650 && y(1) >= 800 && y(1) < 885) && ...
        (x(2) == 650 && y(2) >= 800 && y(2) < 885)) || ...
       ((x(1) > 650 && x(1) <= 800 && y(1) == 800) && ...
        (x(2) > 650 && x(2) <= 800 && y(2) == 800))
        [loc_ER_cal,loc_Cyto_cal] = loc_cal_func(inter_bound_p,loc2glb,ER_cal,Cyto_cal,compartment);
        
        Jserca = Vp*(thap*blocker*CPA_on*loc_Cyto_cal^2-K_bar_serca*loc_ER_cal^2)/...
                    (loc_Cyto_cal^2+k_serca^2);
        if compartment == "ER_cal"
           z = Jserca*buff_scale*gamma;
        elseif compartment == "Cyto_cal"
           z = -1*Jserca*buff_scale;
        else
           fprintf("Something is wrong in gN function (SERCA)");
        end
    % Deep IP3R and SERCA fluxes  
    elseif ((x(1) >= 550 && x(1) <= 1850 && y(1) == 200) && ...
            (x(2) >= 550 && x(2) <= 1850 && y(2) == 200)) || ...
           ((x(1) >= 550 && x(1) <= 1850 && y(1) == 350) && ...
            (x(2) >= 550 && x(2) <= 1850 && y(2) == 350)) || ...
           ((x(1) == 1850 && y(1) >= 200 && y(1) <= 350) && ...
            (x(2) == 1850 && y(2) >= 200 && y(2) <= 350))
            
        deep_Jserca_ratio = 1;
        deep_Jipr_ratio = 1;
        
        [loc_ER_cal,loc_Cyto_cal] = loc_cal_func(inter_bound_p,loc2glb,ER_cal,Cyto_cal,compartment);
        
        Jserca = deep_Jserca_ratio*Vp*(thap*CPA_on*loc_Cyto_cal^2-K_bar_serca*loc_ER_cal^2)/...
                    (loc_Cyto_cal^2+k_serca^2);
        
        loc_h_pointer_1 = (deep_h_p(1,:) == x(1) & deep_h_p(2,:) == y(1));
        loc_h_pointer_2 = (deep_h_p(1,:) == x(2) & deep_h_p(2,:) == y(2));
        h = deep_h;
        
        loc_h_1 = h(loc_h_pointer_1);
        loc_h_2 = h(loc_h_pointer_2);
        
        loc_h_c = (loc_h_1+loc_h_2)/2;
        
        m_alpha = loc_Cyto_cal^4/(Kc^4+loc_Cyto_cal^4);
        m_beta = m_alpha;
        h_bar = Kh^4/(Kh^4+loc_Cyto_cal^4);
        B = Ago^2/(Kp^2+Ago^2);
        A = 1 - B;

        alpha = A*(1-m_alpha*h_bar);
        beta = B*m_beta*loc_h_c;
        P0 = beta/(beta+0.4*(beta+alpha));
        Jipr = deep_Jipr_ratio*kf*P0*(loc_ER_cal-loc_Cyto_cal);
        
        if compartment == "ER_cal"
            z = (Jserca-1*Jipr)*buff_scale*gamma;
        elseif compartment == "Cyto_cal"
            z = (Jipr-Jserca)*buff_scale;
        else
            fprintf("Something is wrong in gN function (IPR)");
        end
    % Cortical IP3R flux
    elseif ((x(1) >= 1300 && x(1) <= 1700 && y(1) == 800) && ...
            (x(2) >= 1300 && x(2) <= 1700 && y(2) == 800))
        [loc_ER_cal,loc_Cyto_cal] = loc_cal_func(inter_bound_p,loc2glb,ER_cal,Cyto_cal,compartment);
    
        loc_h_pointer_1 = (cortical_h_p(1,:) == x(1) & cortical_h_p(2,:) == y(1));
        loc_h_pointer_2 = (cortical_h_p(1,:) == x(2) & cortical_h_p(2,:) == y(2));
        h = cortical_h;
            
        loc_h_1 = h(loc_h_pointer_1);
        loc_h_2 = h(loc_h_pointer_2);
            
        loc_h_c = (loc_h_1+loc_h_2)/2;
            
        m_alpha = loc_Cyto_cal^4/(Kc^4+loc_Cyto_cal^4);
        m_beta = m_alpha;
        h_bar = Kh^4/(Kh^4+loc_Cyto_cal^4);
        B = Ago^2/(Kp^2+Ago^2);
        A = 1 - B;

        alpha = A*(1-m_alpha*h_bar);
        beta = B*m_beta*loc_h_c;
        P0 = beta/(beta+0.4*(beta+alpha));
        Jipr = kf*P0*(loc_ER_cal-loc_Cyto_cal);
            
        if compartment == "ER_cal"
            z = -1*Jipr*buff_scale*gamma;
        elseif compartment == "Cyto_cal"
            z = Jipr*buff_scale;
        else
            fprintf("Something is wrong in gN function (IPR)");
        end
    % SOCE and PMCA fluxes    
    elseif ((x(1) >= 0 && x(1) <= 2350 && y(1) == 900) && ...
            (x(2) >= 0 && x(2) <= 2350 && y(2) == 900))
        if compartment == "Cyto_cal"
            loc_Cyto_cal_1 = Cyto_cal(loc2glb(1,:));
            loc_Cyto_cal_2 = Cyto_cal(loc2glb(2,:));
        else
            fprintf("Something is wrong with Jpm");
        end
            
        loc_Cyto_cal = (loc_Cyto_cal_1+loc_Cyto_cal_2)/2;
        
        Jpm = Vpm*loc_Cyto_cal^2/(Kpm^2+loc_Cyto_cal^2);
        
        % SOCE location
        if ((x(1) >= 550) && (x(1) <= 650) && (y(1) == 900) && ...
            (x(2) >= 550) && (x(2) <= 650) && (y(2) == 900))
                
            STIM = K_stim^STIM_Hill/(K_stim^STIM_Hill+mean(ER_cal,"all")^STIM_Hill);

            Jin = alpha1*(STIM);
                
            if compartment == "Cyto_cal"
                z = (ext_cal_on*(Jin+alpha0)-1*Jpm)*buff_scale;
            else
                fprintf("Something is wrong in gN function (Jin)");
            end
        else
            if compartment == "Cyto_cal"
                z = (alpha0*ext_cal_on-1*Jpm)*buff_scale;
            else
                fprintf("Something is wrong in gN function (Jpm)");
            end
        end
    else
        z = 0;
    end
end

% Function to calculate h variable at each time step
function h = calculate_h_new(h_p,Cyto_cal,h_old,t_step)
    Kh = 0.8;
    Ktau = 0.2;
    Tau_max = 200;
    
    h = zeros(size(h_p,2),1);
    %loc_Cyto_cal = zeros(size(h_p,2));
    for i = 1:size(h_p,2)
        loc_Cyto_cal = Cyto_cal(h_p(3,i));

        h_inf = Kh^4/(Kh^4+loc_Cyto_cal^4);
        Tau_h = Tau_max*Ktau^4/(Ktau^4+loc_Cyto_cal^4);

        h_hold = t_step*(h_inf - h_old(i))/Tau_h + h_old(i);
        if h_hold < 1e-27
            h_hold = 0;
        elseif h_hold > 1
            h_hold = 1;
        end
        h(i) = h_hold;
    end
end

%% Functions for FEM
% These functions follow Larson and Bengzon 2013 book with slight
% modifications to account for interior boundary condition

% A useful function that helps to calculate the Stiffness matrix
function [area,b,c] = HatGradients(x,y)
    area = polyarea(x,y);
    b = [y(2)-y(3);
        y(3)-y(1);
        y(1)-y(2)]/2/area;
    c = [x(3)-x(2);
        x(1)-x(3);
        x(2)-x(1)]/2/area;
end

%Mass matrix function
function M = MassAssembler2D(p,t)
    np = size(p,2); %set the number of nodal points into np
    nt = size(t,2); %set the number of triangle mesh into nt
    M = sparse(np,np); %allocate mass matrix
    for K = 1:nt
        % local to global map
        loc2glb = t(1:3,K); % Getting the nodal info from row 5 to 7 (for subdomain) instead of 1 to 3 (original index)
        x = p(1,loc2glb); % node x coordinates
        y = p(2,loc2glb); % node y coordinates
        area = polyarea(x,y);
        MK = [2 1 1;
            1 2 1;
            1 1 2]/12*area;
        M(loc2glb,loc2glb) = M(loc2glb,loc2glb) + MK;
    end
end

%Stiffness matrix function
function A = StiffnessAssembler2D(p,t,Diff,compartment)
    np = size(p,2); %set the number of points into np
    nt = size(t,2); %set the number of triangle mesh into nt
    A = sparse(np,np); %allocate stiffness matrix
    for K = 1:nt
        %local to global map
        loc2glb = t(1:3,K); % Getting the nodal info from row 5 to 7 (for subdomain) instead of 1 to 3 (original index)
        x = p(1,loc2glb); % node x coordinates
        y = p(2,loc2glb); % node y coordinates
        [area,b,c] = HatGradients(x,y);
        avg_Diff = mean(Diff(x,y,compartment));
        AK = avg_Diff*(b*b'+c*c')*area; %stiffness matrix's element
        %add local stiffness element to global stiffness
        A(loc2glb,loc2glb) = A(loc2glb,loc2glb) + AK;
    end
end

%Additional matrix to be added to the Stiffness matrix with the given
%boundary conditions
function R = RobinStiff2D(p,e,kappa,compartment)
    np = size(p,2); %number of nodes
    ne = size(e,2); %number of boundary edges
    R = sparse(np,np); %allocate boundary matrix
    for E = 1:ne
        %boundary nodes (local to global mapping)
        loc2glb = e(1:2,E); %Getting the nodal info from row 8 and 9 (for subdomain) instead of row 1 and 2
        x = p(1,loc2glb);
        y = p(2,loc2glb);
        len = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2); %edge length
        xc = mean(x);
        yc = mean(y);
        k = kappa(xc,yc,compartment);
        RE = k/6*[2 1; 1 2]*len; %edge boundary matrix
        R(loc2glb,loc2glb) = R(loc2glb,loc2glb) + RE;
    end
end

%Additional vector to be added to the Load vector with the given
%boundary conditions
function r = RobinLoad2D(p,e,kappa,gD,gN,...
    ER_cal,Cyto_cal,inter_bound_p,compartment, ...
    cortical_h_p,cortical_h,deep_h_p,deep_h,i,t_step)
    
    np = size(p,2); %number of nodes
    ne = size(e,2); %number of boundary edges
    r = zeros(np,1);

    for E = 1:ne
        loc2glb = e(1:2,E); %Getting the nodal info from row 8 and 9 (for subdomain) instead of row 1 and 2
        x = p(1,loc2glb);
        y = p(2,loc2glb);
        len = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2); %edge length
        xc = mean(x);
        yc = mean(y);
        
        res_gN = gN(x,y,loc2glb,...
            ER_cal,Cyto_cal,inter_bound_p,compartment, ...
            cortical_h_p,cortical_h,deep_h_p,deep_h,i,t_step);
        tmp = kappa(xc,yc,compartment)*gD(xc,yc,compartment)+res_gN;
        
        rE = tmp*[1; 1]*len/2;
        r(loc2glb) = r(loc2glb) + rE;
    end
end

% The Load vector function
function L = LoadAssembler2D(p,t,f,compartment)
    np = size(p,2);
    nt = size(t,2);
    L = zeros(np,1);

    for K = 1:nt
        %local to global map
        loc2glb = t(1:3,K); % Getting the nodal info from row 5 to 7 (for subdomain) instead of 1 to 3 (original index)
        x = p(1,loc2glb); % node x coordinates
        y = p(2,loc2glb); % node y coordinates
        area = polyarea(x,y);
        LK = [f(x(1),y(1),compartment);
            f(x(2),y(2),compartment);
            f(x(3),y(3),compartment)]/3*area; %local load vector element based on corner quadrature formula
        %add local load element to global load
        L(loc2glb) = L(loc2glb) + LK;
    end
end

