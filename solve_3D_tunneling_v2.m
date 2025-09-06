

% Test code to solve a Ca model in a 3D domain with separated ER and
% cytoplasm

% Solve Ca diffusion model, with an internal ER compartment, 
% SERCA flux and IPR on the ER membrane, and Jin and Jout on the PM. 

% This version assumes that the SOCE flux goes directly into the ER. So
% it's a simplified version of the model with a really large SERCA flux
% close to the SOCE channels.

% This version also simplifies the ce-c (in JIPR) to just ce. Will be a
% very good approximation in the physiological regime. This change helps
% with the scaling of ce.

% YOU HAVE TO RUN PREPARE_MESH.M FIRST, TO GET THE MESH READY.

clear all
close all
clc

global s_PM s_ER v_cyt v_ER 
global np_c np_e np_h cyt_node cyt_node_PM
global PMtriarea PMtriarea_SOC ERtriarea_IPR ERtriarea_SERCA 
global tets_cyt_hom tets_ER_hom
global tetvol_cyt tetvol_ER

global tt
global par


load('mesh.mat')
getparams

f1 = figure(1);
f1.Position=[1000,1000,500,500];
f2 = figure(2);
f2.Position=[1500,1000,500,500];

np_c = size(v_cyt.POS,1);
np_e = size(v_ER.POS,1);
np_h = np_c;                        % Cheating. Solving h at every c node, which is very inefficient. But the c and h nodes always match.


%% Initial conditions
% Initial condition. Note the structure.
% First np_c entries are c
% Next np_e entries are e
% Next np_h entries are h. So h is defined only on the ER surface nodes.

if par.initial_choose == 0
     xnew(1:np_c) = 0.001;
     xnew(np_c+1:np_c+np_e) = 10;
     xnew(np_c+np_e+1:np_c+np_e+np_h) = 0.99;
     xnew = xnew';
elseif par.initial_choose == 1
    load initial_rest.mat
elseif par.initial_choose == 2
    load initial_short_recovery.mat
elseif par.initial_choose == 3
     xnew(1:np_c) = 0.03;
     xnew(np_c+1:np_c+np_e) = 300;
     xnew(np_c+np_e+1:np_c+np_e+np_h) = 0.6;
     xnew = xnew';
elseif par.initial_choose == 4
    load initial.mat
end

keepc(1,:) = xnew(1:np_c);
keepe(1,:) = xnew(np_c+1:np_c+np_e);



%% Mass and stiffness matrices
% first run prepare_mesh to load the triangles, tets, and precompute areas
% etc. Load up the result here.

ntets_c = size(v_cyt.TETS,1);
ntets_e = size(v_ER.TETS,1);

[mass_c,stiff_c] = makematrices(v_cyt.POS,np_c,v_cyt.TETS,ntets_c,tetvol_cyt,par.Dc);            % Get the mass and stiffness matrices for the cytoplasm                                   
[mass_e,stiff_e] = makematrices(v_ER.POS,np_e,v_ER.TETS,ntets_e,tetvol_ER,par.De);               % Get the mass and stiffness matrices for the ER
% the combined matrices, to solve for c and ce together
mass = [mass_c       zeros(np_c,np_e);
       zeros(np_e,np_c)   mass_e];
stiff = [stiff_c       zeros(np_c,np_e);
       zeros(np_e,np_c)   stiff_e];

%% the solve

time(1) = 0;
par.time = time(1);
tic
% Now do the actual solve. The solves for the cytoplasm and the ER must be
% done at the same time, as the ER boundary fluxes involve terms in both c
% and ce.

i=1;
while time(i)<par.tend
    par.time = time(i);
    tt = time(i);               % for passing the time to all the fluxes, to get time-dependent behaviour
    if (min(xnew<0) || max(xnew>5000))
        fprintf('Crap here, at %4.2f\n',tt)
        return
    end
    %fprintf('max c is %4.2f\n',max(xnew(1:np_c)))
    %fprintf('min c is %4.2f\n\n',min(xnew(1:np_c)))
    xold = xnew;
    react = get_load(xold);                         % make the load vector

    % determine the time step. This is a cheap and nasty adaptive time step
    % method. It could be done a lot lot better. Without this adaptive time
    % stepping, the computation doesn't really work.
    if time(i)<0.05
        delt = 0.001;   % force an initial small time step
    else
        delt = par.adapt_factor*max(abs(xnew(1:np_c)))/max(abs(react));
    end

    time(i+1) = time(i) + delt;
                   % easiest way to pass time to the subroutines
    fprintf('%5.4f completed; sample ER value is %5.1f\n',time(i)/par.tend,xnew(np_c+1))                                    % For tracking how the computation is going
    Amat_c = (mass_c + delt*stiff_c);        
    Amat_e = (mass_e + delt*stiff_e);
    Amat = [Amat_c       zeros(np_c,np_e);
       zeros(np_e,np_c)   Amat_e];                  % Amat contains only entries for the diffusing variables

    % do the time step
    xnew(1:np_c+np_e) = Amat\(mass*xold(1:np_c+np_e)+ delt*react(1:np_c+np_e));          % c and ce are solved by matrix inversion, as they diffuse.
    xnew(np_c+np_e+1:np_c+np_e+np_h) = xold(np_c+np_e+1:np_c+np_e+np_h) + delt*react(np_c+np_e+1:np_c+np_e+np_h); % Forward Euler for h.

    % These are for if you need to do the total calcium check
    %total_c = computetotalcalcium(xnew(1:np_c),np_c,tets_c,ntets_c,tetvol_c);
    %total_e = computetotalcalcium(xnew(np_c+1:np_c+np_e),np_e,tets_e,ntets_e,tetvol_e);
    %total(i) = total_c + total_e;

    keepc(i+1,:) = xnew(1:np_c);
    keepe(i+1,:) = xnew(np_c+1:np_c+np_e);
    delt_keep(i) = delt;                        % Just for curiosity, to see how the time step is changing.

    i = i+1;

    %figure(1)
    %plot(time,keepc,'LineWidth',2)
    %figure(2)
    %plot(time,keepe(:,ERIsp),'LineWidth',2)
end
toc

save output_temp1.mat
%save initial.mat xnew

%% plots
figure(1)
plot(time,keepc)
%ylim([0,0.5])
figure(2)
plot(time,keepe)
figure(3)
plot(delt_keep)
%figure(3)
%plot(total)



%% ----------------------------------------------------------------------------

%% check by computing total calcium

function out = computetotalcalcium(x,np,tets,ntets,tetvol)

total = 0;
for i=1:ntets
    loc2glb = tets(i,1:4); % local-to-global map
    av = (x(loc2glb(1)) + x(loc2glb(2)) + x(loc2glb(3)) + x(loc2glb(4)) )/4 * tetvol(i);
    total = total + av;
end

out = total;
end

%% mass and stiffness matrices
function [mass,stiff] = makematrices(p,np,tets,ntets,tetvol,D)
 
    stiffc=sparse(np,np);
    small_mass = sparse(np,np); % allocate mass matrix
    
    for K = 1:ntets % loop over elements
        loc2glb = tets(K,1:4); % local-to-global map
        P1 = p(loc2glb(1),:); P2 = p(loc2glb(2),:); P3 = p(loc2glb(3),:);   P4 = p(loc2glb(4),:);   % mostly just for convenience
        TT = [1 P1; 1 P2; 1 P3; 1 P4];
        %volume = abs(det(TT))/6;
        
        % first the mass matrix
        MK = [2 1 1 1;
        1 2 1 1;
        1 1 2 1;
        1 1 1 2]/20*tetvol(K); % element mass matrix
        small_mass(loc2glb,loc2glb) = small_mass(loc2glb,loc2glb) + MK; % add element masses to M
        
        % then the stiffness matrix
        dum = inv(TT'); b = dum(:,2); c = dum(:,3); d = dum(:,4);   % This method of getting the coeffs follows Gockenbach
        Ac = D*(b*b' + c*c' + d*d')*tetvol(K);                      % element stiffness matrix
        stiffc(loc2glb,loc2glb) = stiffc(loc2glb,loc2glb) + Ac;     % add element stiffnesses to global stiffness matrix
    end
    
    % Since h doesn't diffuse, it is solved directly, and
    % doesn't need entries in the mass and stiffness matrices.
    mass  = sparse(small_mass);
    stiff = sparse(stiffc);

end

%% load vector
function out=get_load(u)
    global PMtriarea PMtriarea_SOC ERtriarea_IPR ERtriarea_SERCA 
    global s_PM s_ER v_cyt v_ER 
    global par np_c np_e np_h cyt_node cyt_node_PM
    global tets_cyt_hom tets_ER_hom
    global tetvol_cyt tetvol_ER

    %   First np_c components of u    -   c
    %   Next np_e components of u    -   e
    %   Next np_h components of u    -   h
    
    c = u(1:np_c);
    e = u(np_c+1:np_c+np_e);
    h = u(np_c+np_e+1:np_c+np_e+np_h);

    eav = mean(e);
    
    load_c = zeros(np_c,1);
    load_e = zeros(np_e,1);
    
    % -------------------------------------

    % We now have to include all the reaction terms, but ONLY in the
    % homogenised region (so that each node has both c and ce defined
    % there). There is no explicit boundary surface  between the homogenised region
    % and the microdomain; the difference is handled simply by including
    % the reaction terms in the homogenised region only, and (if desired)
    % by changing the diffusion coefficients in the homogenised region.


%   v_cyt.TETS(tets_cyt_hom(i),1:4) is the same tet as
%   v_ER.TETS(tets_ER_hom(i),1:4)

    nhtets = size(v_cyt.TETS(tets_cyt_hom,1:4),1);
    for K = 1:nhtets
        loc2glb_c = v_cyt.TETS(tets_cyt_hom(K),1:4);
        loc2glb_e = v_ER.TETS(tets_ER_hom(K),1:4);

        % get the reactions at each of the homogenised region tetrahedron nodes. 
        r1 = gethomogreactions(c(loc2glb_c(1)),e(loc2glb_e(1)),h(loc2glb_c(1)),par);
        r2 = gethomogreactions(c(loc2glb_c(2)),e(loc2glb_e(2)),h(loc2glb_c(2)),par);
        r3 = gethomogreactions(c(loc2glb_c(3)),e(loc2glb_e(3)),h(loc2glb_c(3)),par);
        r4 = gethomogreactions(c(loc2glb_c(4)),e(loc2glb_e(4)),h(loc2glb_c(4)),par);

        b_cal = [r1(1); r2(1); r3(1); r4(1)]/4*tetvol_cyt(tets_cyt_hom(K));         % local element load vector.
        load_c(loc2glb_c) = load_c(loc2glb_c) + b_cal;                              % add element loads to global load vector
        load_e(loc2glb_e) = load_e(loc2glb_e) - par.gamma*b_cal;                    % add element loads to global load vector
    end
     
    % -------------------------------------
    
    % Now add in the boundary conditions. These are obtained by integrating the
    % appropriate flux over the surface triangles, and adding these integrals to the
    % appropriate entries in the load vector.
    
    % There are two lots of boundary conditions; one at the PM, the other at
    % the ER membrane. Also, each boundary condition is split into two
    % bits, as the SERCA, IPR, SOC and PM pump fluxes all occur on
    % different boundary triangles.

    % Note that you have to be very careful about the node where each flux
    % occurs, because the boundary nodes are labelled differently from the
    % cyt nodes. Hence the use of cyt_node_PM.

    % Node i in the s_PM list corresponds to node cyt_node_PM(i) in the cyt
    % list.
    
    % First do the PM pump flux. This occurs on every PMCA triangle
    for K = 1:size(s_PM.TRIANGLES_PMCA,1)
        loc2glb = s_PM.TRIANGLES_PMCA(K,1:3);
        r1 = getPMbndyreact_pump(c(cyt_node_PM(loc2glb(1))),par);
        r2 = getPMbndyreact_pump(c(cyt_node_PM(loc2glb(2))),par);
        r3 = getPMbndyreact_pump(c(cyt_node_PM(loc2glb(3))),par);
        
        bndy = [r1; r2; r3]/3*PMtriarea(K);
        load_c(cyt_node_PM(loc2glb)) = load_c(cyt_node_PM(loc2glb)) + bndy;
    end

    % Next do the SOCE flux. This occurs only on the SOCE triangles (on the
    % PM). It also depends only on the average ER Ca concentration, which
    % is a total cheat, but doing anything else would be far too difficult.

    % In this version of the model there is no SOCE flux into the
    % cytoplasm. It is assumed to go directly into the cortical ER
    % instead. 

            % for K = 1:size(s_PM.TRIANGLES_SOC,1)
            %     loc2glb = s_PM.TRIANGLES_SOC(K,1:3);
            %     r1 = getPMbndyreact_soce(eav,par);
            %     r2 = getPMbndyreact_soce(eav,par);
            %     r3 = getPMbndyreact_soce(eav,par);
            % 
            %     bndy = [r1; r2; r3]/3*PMtriarea_SOC(K);
            %     load_c(loc2glb) = load_c(loc2glb) + bndy;
            % end

    
    % The ER surface must be done more carefully, as it separates two
    % regions, each of which has a variable for which we are solving. This
    % means that every triangle on the ER is in two different tets; one in
    % the ER, one in the cytoplasm. But these tets have DIFFERENT
    % numberings in general. So, for each node on the ER surface, you need to know
    % its number in the ER list AND its number in the cytoplasm list. This correspondence is
    % calculated in prepare_mesh.m. 
    
    % Node i in the ER surface list corresponds to node cyt_node(i) in the cytoplasm list.

    % Furthermore, the boundary fluxes occur in two different places (SERCA
    % and IPR) so they need to be added in separately.

    % Note that here we are assuming the SAME buffering capacity in the ER
    % and the cytoplasm, which is probably inaccurate. One of these
    % fluxes should probably be multiplied by some factor. Not, maybe, as
    % big as par.gamma, which takes care of this problem in the homogenised
    % region, but it should probably still be there.

    for K = 1:size(s_ER.TRIANGLES_IPR,1)
        loc2glb = s_ER.TRIANGLES_IPR(K,1:3);   % these are the nodes in the ER list
        r1 = get_react_ipr(c(cyt_node(loc2glb(1))),e(loc2glb(1)),h(cyt_node(loc2glb(1))),par,par.k_f);
        r2 = get_react_ipr(c(cyt_node(loc2glb(2))),e(loc2glb(2)),h(cyt_node(loc2glb(2))),par,par.k_f);
        r3 = get_react_ipr(c(cyt_node(loc2glb(3))),e(loc2glb(3)),h(cyt_node(loc2glb(3))),par,par.k_f);
        
        bndy = [r1; r2; r3]/3*ERtriarea_IPR(K);
        load_e(loc2glb) = load_e(loc2glb) + bndy;
        load_c(cyt_node(loc2glb)) = load_c(cyt_node(loc2glb)) - bndy;
    end

    for K = 1:size(s_ER.TRIANGLES_SERCA,1)
        loc2glb = s_ER.TRIANGLES_SERCA(K,1:3);   % these are the nodes in the ER list
        r1 = get_react_serca(c(cyt_node(loc2glb(1))),e(loc2glb(1)),h(cyt_node(loc2glb(1))),par,par.Vs);
        r2 = get_react_serca(c(cyt_node(loc2glb(2))),e(loc2glb(2)),h(cyt_node(loc2glb(2))),par,par.Vs);
        r3 = get_react_serca(c(cyt_node(loc2glb(3))),e(loc2glb(3)),h(cyt_node(loc2glb(3))),par,par.Vs);
        
        bndy = [r1; r2; r3]/3*ERtriarea_SERCA(K);
        load_e(loc2glb) = load_e(loc2glb) + bndy;
        load_c(cyt_node(loc2glb)) = load_c(cyt_node(loc2glb)) - bndy;
    end

    % Now add the additional SOCE flux for these ER triangles, but add it
    % ONLY to the ER flux, not to the cytoplasm flux.

    for K = 1:size(s_ER.TRIANGLES_SERCA,1)
        loc2glb = s_ER.TRIANGLES_SERCA(K,1:3);   % these are the nodes in the ER list
        r1 = getPMbndyreact_soce(eav,par);
        r2 = getPMbndyreact_soce(eav,par);
        r3 = getPMbndyreact_soce(eav,par);
        
        bndy = [r1; r2; r3]/3*ERtriarea_SERCA(K);
        load_e(loc2glb) = load_e(loc2glb) + bndy;
    end
    
    out(1:np_c+np_e,1)             =    [load_c; load_e];       
     
    % The variables with no diffusion can be done differently. No loop over
    % the triangles is needed, just a simpler loop over the nodes, but then
    % time must be stepped forward using a forward Euler, not with a solve
    % using the mass matrix.
     
    for K=1:np_h
        out(np_c+np_e+K,1)   =  get_h_react(c(K),h(K),par);      % Note how this is done for ALL the c nodes. Rather inefficient.      
    end

end


%% ODE for h
function out = get_h_react(c,h,par)  

    % IPR
    H_inf = par.K_h^4 / ( par.K_h^4 + c^4 );
    TAU = par.tau_max*par.K_tau^4/(par.K_tau^4+c^4);
    out = (H_inf - h) / TAU;

end

%% PM pump fluxes on PM boundary
function out = getPMbndyreact_pump(c,par)  
    Jpm = par.Vpm*c*c/(par.Kpm*par.Kpm + c*c);
    out = -par.delta*Jpm;                           % Note the sign of out.
end

%% SOC fluxes on PM boundary
function out = getPMbndyreact_soce(eav,par)  
    if par.ext_Ca
        Jsoce = par.alpha0 + par.Vsoc*(par.Ksoc^4)/(eav^4 + par.Ksoc^4);
    else
        Jsoce = 0; % 0 external Ca case
    end
    out = par.delta*Jsoce;
end

%% IPR fluxes on ER boundary
function out = get_react_ipr(c,e,h,par,IPRdensity)

    if (par.time<par.t_stim)
        par.ip = 0;
    else
        par.ip = par.ip_stim;
    end

    phi_c = c^4 / ( c^4 + par.K_c^4 );
    phi_p = par.ip^2 / ( par.K_p^2 + par.ip^2);
    phi_p_down = par.K_p^2 / ( par.K_p^2 + par.ip^2);
    H_inf = par.K_h^4 / ( par.K_h^4 + c^4 );
    TAU = par.tau_max*par.K_tau^4/(par.K_tau^4+c^4);
    beta = phi_p * phi_c * h;
    alpha = phi_p_down * ( 1 - phi_c * H_inf );
    po = beta/(beta+par.k_beta*(beta+alpha));

    JIPR = IPRdensity*po*(e);                       % simplified version
    out = -JIPR;                                    % Note the sign of out.
        
end
%% SERCA fluxes on ER boundary 
function out = get_react_serca(c,e,~,par,serca_density)  
    
    if par.CPA
        Jserca = serca_density*(- par.Kbar*e*e)/(par.Ks*par.Ks + c*c);   % With CPA
    else
        Jserca = serca_density*(c*c - par.Kbar*e*e)/(par.Ks*par.Ks + c*c); % No CPA
    end
    out = Jserca;
        
end

%% Reactions in the homogenised region
function out = gethomogreactions(c,e,h,par)  

    Jserca = get_react_serca(c,e,h,par,par.Vs_homog);
    JIPR = -get_react_ipr(c,e,h,par,par.k_f_homog);          % The minus sign because of the minus sign in getERbndyreact_ipr
    out = JIPR - Jserca;
        
end

%% other parameters
function getparams
global par

    % 0 for empty. (manual)
    % 1 for rest. (initial_rest.mat) 
    % 2 for short recovery. 
    % 3 for manual rest.
    % 4 for initial.mat
    par.initial_choose = 0; 

    par.ip_stim = 1;
    par.t_stim = 0.25;
    par.tend = 2;
        
    par.adapt_factor = 0.01;

    par.ext_Ca = true;          % true for external Ca, false for no external Ca 
    par.CPA = false;            % true for CPA, false for no CPA
    par.lam = 60;               % For rescaling the ER

    % ER fluxes
    par.k_f = 20/par.lam;
    par.k_f_homog = 0.1/par.lam;

    % Don't change stuff below here
    
    par.delta = 5; 

    par.Dc = 5;
    par.De = 5;

    par.Vs = 10;
    par.Vs_homog = 0.5;

    par.Ks = 0.2;
    par.Kbar = 0.00001957/(par.lam*par.lam);
    par.gamma = 8*par.lam;

    % PM fluxes
    par.Vpm = 0.1;    
    par.Vsoc = 20;
    par.alpha0 = 0.0001;
    par.Ksoc = 10*par.lam;
    par.Kpm = 0.3;
         
    % IPR parameters
    par.K_c = 0.2;
    par.K_p = 0.2;
    par.K_h = 0.08;
    par.tau_max = 10;
    par.K_tau = 0.1;
    par.k_beta = 0.4;

end