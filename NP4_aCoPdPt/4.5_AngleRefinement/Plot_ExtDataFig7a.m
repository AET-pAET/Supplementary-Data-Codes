clear all; close all
%%
angles_GT=load('CoPdPt_Angle_Dose1p7e4.mat');% Ground truth
angles_GT=angles_GT.angles_refined_out;
angles      = importdata('CoPdPt_ExperimentalAngles.mat');% Angle from experiment
angles_OPT      = importdata('CoPdPt_RefinedAngles.mat'); % Refined Angles
angles_OPT=angles_OPT.Angle_OPT;
%%
% compute rotation matrix R w.r.t euler angles {phi,theta,psi}
% rotmat1 = MatrixQuaternionRot(vec1,phi);
% rotmat2 = MatrixQuaternionRot(vec2,theta);
% rotmat3 = MatrixQuaternionRot(vec3,psi);

colors = [0 0 160
    0 102 102
    220 111 0
    237 67 62]/255;

ImagNo=[1:55];

figure; clf; set(gcf,'Position',[0,0,750,190]);
subplot(1,3,1);
plot(angles(:,1)-angles_GT(:,1),'+','markersize',6,'color',colors(3,:),'MarkerFaceColor',[1 1 1],'LineWidth',1.0)
hold on
plot(angles_OPT(:,1)-angles_GT(:,1),'.','markersize',6,'color',colors(1,:),'MarkerFaceColor', [1, 1, 1],'LineWidth',1.0)%  [0.3, 0.3, 0.8]
box on
set(gca,'Layer','top')
h = findobj(gca,'Type','patch');
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1);
xlabel('Image Number');
ylabel('Angle error (°)');
lgd = legend('Before refinement','After refinement');% 
lgd.Box = 'off'; 
set(legend,'FontSize',8)
lgd.Orientation = "horizontal";
lgd.NumColumns = 1;
lgd.Location= 'northwest'
ylim([-4 5])
lgd.ItemTokenSize = [20, 20];


subplot(1,3,2);
plot(angles(:,2)-angles_GT(:,2),'+','markersize',6,'color',colors(3,:),'MarkerFaceColor',[1 1 1],'LineWidth',1.0)
hold on
plot(angles_OPT(:,2)-angles_GT(:,2),'.','markersize',6,'color',colors(1,:),'MarkerFaceColor', [1, 1, 1],'LineWidth',1.0)% 
box on
set(gca,'Layer','top')
h = findobj(gca,'Type','patch');
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1);
xlabel('Image No.');
ylabel('Angle error (°)');
lgd = legend('Before refinement','After refinement');% 
lgd.Box = 'off'; 
set(legend,'FontSize',8)
lgd.Orientation = "horizontal";
lgd.NumColumns = 1;
lgd.Location= 'northwest'
ylim([-4 5])
lgd.ItemTokenSize = [20, 20];


subplot(1,3,3);
plot(angles(:,3)-angles_GT(:,3),'+','markersize',6,'color',colors(3,:),'MarkerFaceColor',[1 1 1],'LineWidth',1.0)
hold on
plot(angles_OPT(:,3)-angles_GT(:,3),'.','markersize',6,'color',colors(1,:),'MarkerFaceColor', [1, 1, 1],'LineWidth',1.0)% 
box on
set(gca,'Layer','top')
h = findobj(gca,'Type','patch');
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1);
xlabel('Image No.');
ylabel('Angle error (°)');
lgd = legend('Before refinement','After refinement');% 
lgd.Box = 'off'; 
set(legend,'FontSize',8)
lgd.Orientation = "horizontal";
lgd.NumColumns = 1;
lgd.Location= 'northwest'
ylim([-4 5])
lgd.ItemTokenSize = [20, 20];
