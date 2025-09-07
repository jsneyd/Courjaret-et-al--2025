
% For having a look at the output from solve_3D runs. Choose the run you
% want to inspect.


clear all
clc
close all

load output_tunneling_Dc5_De5.mat


eav = mean(keepe(:,ERIsp),2);

%%
figure(1);
trisurf(s_ER.TRIANGLES(:,1:3),s_ER.POS(:,1),s_ER.POS(:,2),s_ER.POS(:,3),'FaceColor','cyan','FaceAlpha',0.5,'LineStyle',"none");
hold on
trisurf(s_PM.TRIANGLES(:,1:3),s_PM.POS(:,1),s_PM.POS(:,2),s_PM.POS(:,3),'FaceColor','red','FaceAlpha',0.1,'LineStyle',"none");
axis equal
% trisurf(s_PM.TRIANGLES_SOC,s_PM.POS(:,1),s_PM.POS(:,2),s_PM.POS(:,3),'FaceColor','green','FaceAlpha',0.5,'LineStyle',"none");
% trisurf(s_ER.TRIANGLES_SERCA,s_ER.POS(:,1),s_ER.POS(:,2),s_ER.POS(:,3),'FaceColor','red','FaceAlpha',0.5,'LineStyle',"none");
% trisurf(s_ER.TRIANGLES_IPR,s_ER.POS(:,1),s_ER.POS(:,2),s_ER.POS(:,3),'FaceColor','blue','FaceAlpha',0.5,'LineStyle',"none");

scatter3(v_cyt.POS(Isp,1),v_cyt.POS(Isp,2),v_cyt.POS(Isp,3),75,"filled","red")
text(v_cyt.POS(Isp,1),v_cyt.POS(Isp,2),v_cyt.POS(Isp,3),sprintfc(' %d',1:numel(v_cyt.POS(Isp,1))))
scatter3(v_ER.POS(ERIsp,1),v_ER.POS(ERIsp,2),v_ER.POS(ERIsp,3),75,"filled","black")
text(v_ER.POS(ERIsp,1),v_ER.POS(ERIsp,2),v_ER.POS(ERIsp,3),sprintfc(' %d',1:numel(v_ER.POS(ERIsp,1))))
%set(gca,'XTickLabel',[], 'YTickLabel', [], 'ZTickLabel', [])

[x y] = meshgrid(-1:1); % Generate x and y data
z = -0.3*ones(size(x, 1)); % Generate z data
surf(x, y, z,'FaceAlpha',0.1,'LineStyle',"none") % Plot the surface

ylim([-0.4,0.6])
xlim([-0.8,0.8])
ax = gca;
ax.TickLength = [0 0];

%%
f2=figure(2);
f2.Position=[500,500,1000,350];
subplot(1,2,1)
plot(time,keepc(:,Isp),'LineWidth',2)
legend("cyto 1 (red dots)","cyto 2","cyto 3","cyto 4")
xlabel("time (s)")
ylabel("cytoplasmic Ca conc")

subplot(1,2,2)
plot(time,keepe(:,ERIsp),'LineWidth',2)
legend("ER 1 (green dots)","ER 2","ER 3","ER 4")
xlabel("time (s)")
ylabel("ER Ca conc")

igor_plots = [time' keepc(:,Isp) keepe(:,ERIsp)];
writematrix(igor_plots,'for_igor.dat')
writematrix(igor_plots,'for_khaled.csv')

%%
% f3=figure(3);
% f3.Position=[500,500,1000,350];
% subplot(1,2,1)
% plot(time,get_react_serca(keepc(:,Isp),keepe(:,ERIsp),1,par,par.Vs),'LineWidth',2)
% legend("cyto 1 (red dots)","cyto 2","cyto 3","cyto 4")
% xlabel("time (s)")
% ylabel("SERCA pump flux at selected points")
% 
% subplot(1,2,2)
% plot(time,getPMbndyreact_soce(eav,par),'LineWidth',2)
% xlabel("time (s)")
% ylabel("SOC flux")


%%
figure(4);
trisurf(s_ER.TRIANGLES(:,1:3),s_ER.POS(:,1),s_ER.POS(:,2),s_ER.POS(:,3),'FaceColor','cyan','FaceAlpha',0.5);
hold on
trisurf(s_PM.TRIANGLES(:,1:3),s_PM.POS(:,1),s_PM.POS(:,2),s_PM.POS(:,3),'FaceColor','red','FaceAlpha',0.1);
axis equal
trisurf(s_PM.TRIANGLES_SOC,s_PM.POS(:,1),s_PM.POS(:,2),s_PM.POS(:,3),'FaceColor','green','FaceAlpha',0.5);
trisurf(s_ER.TRIANGLES_SERCA,s_ER.POS(:,1),s_ER.POS(:,2),s_ER.POS(:,3),'FaceColor','red','FaceAlpha',0.5);
trisurf(s_ER.TRIANGLES_IPR,s_ER.POS(:,1),s_ER.POS(:,2),s_ER.POS(:,3),'FaceColor','blue','FaceAlpha',0.5);

scatter3(v_cyt.POS(Isp,1),v_cyt.POS(Isp,2),v_cyt.POS(Isp,3),115,"filled","red")
%text(v_cyt.POS(Isp,1),v_cyt.POS(Isp,2),v_cyt.POS(Isp,3),sprintfc(' %d',1:numel(v_cyt.POS(Isp,1))))
scatter3(v_ER.POS(ERIsp,1),v_ER.POS(ERIsp,2),v_ER.POS(ERIsp,3),115,"filled","black")
%text(v_ER.POS(ERIsp,1),v_ER.POS(ERIsp,2),v_ER.POS(ERIsp,3),sprintfc(' %d',1:numel(v_ER.POS(ERIsp,1))))
set(gca,'XTickLabel',[], 'YTickLabel', [], 'ZTickLabel', [])

[x y] = meshgrid(-1:1); % Generate x and y data
z = -0.3*ones(size(x, 1)); % Generate z data
surf(x, y, z,'FaceAlpha',0.1,'LineStyle',"none") % Plot the surface

ylim([-0.3,0.4])
xlim([-0.5,0.6])
zlim([-0.3,0.4])
ax = gca;
ax.TickLength = [0 0];

%% SERCA fluxes on ER boundary 
function out = get_react_serca(c,e,~,par,serca_density)  
    out = par.gamma*serca_density*(c.*c - par.Kbar*e.*e)./(par.Ks*par.Ks + c.*c);     
end

%% SOC flux
function out = getPMbndyreact_soce(eav,par)  
    Jsoce = par.alpha0 + par.Vsoc*(par.Ksoc^4)./(eav.^4 + par.Ksoc^4); 
    out = par.delta*Jsoce;
end