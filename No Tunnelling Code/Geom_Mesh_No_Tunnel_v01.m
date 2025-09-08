%% Mesh Generation for the No Tunnelling Calcium Model with 2 Subdomains (Cytoplasm and ER)
% This is the code that generates the mesh of the Calcium Model (with
% interior boundary). This mesh generation code follows MATLAB's
% PDE Modeler User Guide.

% Clear everything to ensure an empty environment before running 
clear all
close all
clc

% Set Hmax, Hmin and Hgrad values (parameters that determine mesh size and
% growth rate)
Hmax_par = 70; % Maximum triangular mesh edge
Hmin_par = 3; % Minimum triangular mesh edge
Hgrad_par = 1.75; %Specify how fast the mesh size grow (between 1 and 2)
Hedge_par = 3; %Specify how fine the mesh size is, override Hmin_par (use in cortical region)
Hvertex_par = Hedge_par;

%% Create geometry objects
% Create an ER polygon as the base
ER_geom = [2; 10; 550; 1850; 1850; 650; 650; 1850; 1850; 650; 650; 550;
    200; 200; 350; 350; 770; 770; 800; 800; 885; 885];
% Full cell geom
Full_geom = [2; 4; 0; 2350; 2350; 0;
    0; 0; 900; 900;];
Full_geom = [Full_geom; zeros(length(ER_geom) - length(Full_geom),1)];

% Create a matrix to store the ER polygon
ER_geom_matrix = ER_geom;
% Create name for the ER geometry
ER_geom_name = char('ER');
ER_geom_name = ER_geom_name';
% The formula for the resultant ER geometry
ER_geom_formula = 'ER';
% Create the ER geometry object
[ER_dl,ER_bt] = decsg(ER_geom_matrix,ER_geom_formula,ER_geom_name);
% Create ER PDE object
ER_model = createpde();
% Set the ER geometry into the ER PDE model
geometryFromEdges(ER_model,ER_dl);
% Generate the mesh of the ER geometry
ER_mesh = generateMesh(ER_model,'Hmax',Hmax_par,...
    'Hmin',Hmin_par,'Hedge',{[4 5 8 10],Hedge_par},...
    'Hgrad',Hgrad_par,'GeometricOrder','linear');

% Create a matrix to store the ER polygon and the full cell polygon
geom_mat = [ER_geom, Full_geom];
% Create name for the stored polygon
geom_name = char('ER','Full');
geom_name = geom_name';
% The formula to get the cytoplasm geometry by deducting the ER poylgon
% from the full cell polygon
Cyto_geom = 'Full-ER';
% Create the cytoplasm geometry object
[Cyto_dl,Cyto_bt] = decsg(geom_mat,Cyto_geom,geom_name);
% Create cytoplasm PDE object
Cyto_model = createpde();
% Set the cytoplasm geometry into the cytoplasm PDE model
geometryFromEdges(Cyto_model,Cyto_dl);
% Generate the mesh of the cytoplasm geometry
Cyto_mesh = generateMesh(Cyto_model,'Hmax',Hmax_par,...
    'Hmin',Hmin_par,'Hedge',{[3 4 6 8 10],Hedge_par},...
    'Hgrad',Hgrad_par,'GeometricOrder','linear');

%% Plot the resultant geometry object
% Plot the resultant geometry object to view its edges and subdomains.
% Great to check which edge is the interior boundary.
figure(1)
pdegplot(ER_dl,"EdgeLabels","on","VertexLabels","off")
hold on 
pdegplot(Cyto_dl,"EdgeLabels","off","VertexLabels","off")
xlim([-50 2500])
ylim([-50 900])

figure(2)
pdegplot(ER_dl,"EdgeLabels","off","VertexLabels","off") 
hold on 
pdegplot(Cyto_dl,"EdgeLabels","on","VertexLabels","off")
xlim([-50 2500])
ylim([-50 900])

% Plot the resultant geometry with the generated mesh
figure(3)
pdeplot(ER_model,NodeLabels="off");
hold on
pdeplot(Cyto_model,NodeLabels="off");
xlim([-50 2500])
ylim([-50 900])

%% Extract the Point, Edge, and Connectivity matrix of the resultant geometric object
[ER_p,ER_e,ER_t] = meshToPet(ER_model.Mesh);
[Cyto_p,Cyto_e,Cyto_t] = meshToPet(Cyto_model.Mesh);

% Round up the coordinates to 8 decimal points for easier match of the interior
% shared boundaries/edges
ER_p = round(ER_p,8);
Cyto_p = round(Cyto_p,8);

%% Extracting nodal info from Cytoplasm and ER regions which share the same edge
for i = 1:size(Cyto_p,2)
    Cyto_p(3,i) = i;
end

for i = 1:size(ER_p,2)
    ER_p(3,i) = i;
end

% Extract the nodal info of the interior boundaries from ER
inter_bound_ER_e_pointer = (ER_e(5,:) == 1 | ER_e(5,:) == 2 | ER_e(5,:) == 7 |...
   ER_e(5,:) == 4 | ER_e(5,:) == 10 | ER_e(5,:) == 5 | ER_e(5,:) == 6);
inter_bound_ER_e = ER_e(1:2,inter_bound_ER_e_pointer);
inter_bound_ER_e = reshape(inter_bound_ER_e,1,[]);
inter_bound_ER_e = unique(inter_bound_ER_e);
inter_bound_ER_p = zeros(3,size(inter_bound_ER_e,2));
for i = 1:size(inter_bound_ER_e,2)
    inter_bound_ER_p(:,i) = ER_p(:,inter_bound_ER_e(1,i));
end

% Extract the nodal info of the interior boundaries from cytoplasm
inter_bound_Cyto_e_pointer = (Cyto_e(5,:) == 1 | Cyto_e(5,:) == 7 | Cyto_e(5,:) == 11 |...
    Cyto_e(5,:) == 3 | Cyto_e(5,:) == 4 | Cyto_e(5,:) == 10 | Cyto_e(5,:) == 14);
inter_bound_Cyto_e = Cyto_e(1:2,inter_bound_Cyto_e_pointer);
inter_bound_Cyto_e = reshape(inter_bound_Cyto_e,1,[]);
inter_bound_Cyto_e = unique(inter_bound_Cyto_e);
inter_bound_Cyto_p = zeros(3,size(inter_bound_Cyto_e,2));
for i = 1:size(inter_bound_Cyto_e,2)
    inter_bound_Cyto_p(:,i) = Cyto_p(:,inter_bound_Cyto_e(1,i));
end

% Check to ensure the all nodal info of the interior boundaries between ER and
% cytoplasm match. If not, choose a different Hmax, Hmin, and Hgrad
% parameter values to redefine the mesh in both ER and cytoplasm domains.
inter_bound_p = inter_bound_Cyto_p;
for i = 1:size(inter_bound_Cyto_p,2)
    inter_bound_p_pointer = (inter_bound_ER_p(1,:) == inter_bound_Cyto_p(1,i)) & ...
        (inter_bound_ER_p(2,:) == inter_bound_Cyto_p(2,i));
    inter_bound_p(4,i) = inter_bound_ER_p(3,inter_bound_p_pointer);
end

% Locate the nodal info where the IP3R is placed
n_nodal_cor_IPR_pointer = (Cyto_p(1,:) >= 1250 & Cyto_p(1,:) <= 1650 &...
    Cyto_p(2,:) == 800);
n_nodal_cor_IPR = sum(n_nodal_cor_IPR_pointer);
% Locate the nodal info where the SOCE is placed
n_nodal_cor_SOCE_pointer = (Cyto_p(1,:) >= 500 & Cyto_p(1,:) <= 600 &...
    Cyto_p(2,:) == 900);
n_nodal_cor_SOCE = sum(n_nodal_cor_SOCE_pointer);
% Locate the nodal info where the SERCA pumps are placed
n_nodal_cor_SERCA_pointer = (Cyto_p(1,:) >= 450 & Cyto_p(1,:) <= 900 &...
    Cyto_p(2,:) >= 770 & Cyto_p(2,:) < 885) | (Cyto_p(1,:) >= 450 & Cyto_p(1,:) <= 500 &...
    Cyto_p(2,:) == 885) | (Cyto_p(1,:) >= 550 & Cyto_p(1,:) <= 600 &...
    Cyto_p(2,:) == 885);
n_nodal_cor_SERCA = sum(n_nodal_cor_SERCA_pointer);
% Locate the nodal info where the PMCA is placed
n_nodal_cor_PM_pointer = (Cyto_p(1,:) >= 0 & Cyto_p(1,:) <= 2350 &...
    Cyto_p(2,:) == 900);
n_nodal_cor_PM = sum(n_nodal_cor_PM_pointer);

%% Important info to use
% [ER_p,ER_e,ER_t], [Cyto_p,Cyto_e,Cyto_t], inter_bound_p

save("Mesh_no_tunnel.mat")
