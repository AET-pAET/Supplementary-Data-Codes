close all; clear all;
%%
addpath('1.1_projections');
addpath('1.2_reconstructions');

%%
RESIRE_Dose8p0e3=load('Resire_aSi2nm_FosCetr_Dose8p0e3_cov35pm_BM3D_S0p25N200.mat')
rec_Dose8p0e3=RESIRE_Dose8p0e3.OBJ.reconstruction;

figure
clf; 
set(gcf,'color','w');
set(gcf,'outerposition',[100 100 300 300])
isosurface(rec_Dose8p0e3,0.0112)
axis image;
xlim([64-50,65+50]); ylim([64-50,65+50]); zlim([64-50,65+55]);
xticks(20:20:100); yticks(20:20:100); zticks(20:20:100);
colormap hot
axis tight
axis equal
ax = gca;
ax.LineWidth = 1.0;
set(gca,'linewidth',1.0,'fontsize',11,'fontname','Arial');
[caz,cel] = view()
%print('-r500','-dpdf','si_multislice_defocus_0_8e4_ct_iso_0_3.pdf');
%%
RESIRE_Dose1p6e4=load('Resire_aSi2nm_FosCetr_Dose1p6e4_cov35pm_BM3D_S0p25N200.mat')
rec_Dose1p6e4=RESIRE_Dose1p6e4.OBJ.reconstruction;


figure
clf; 
set(gcf,'color','w');
set(gcf,'outerposition',[100 100 300 300])
isosurface(rec_Dose1p6e4,0.026)
axis image;
xlim([64-50,65+50]); ylim([64-50,65+50]); zlim([64-50,65+55]);
xticks(20:20:100); yticks(20:20:100); zticks(20:20:100);
colormap hot
axis tight
axis equal
ax = gca;
ax.LineWidth = 1.0;
set(gca,'linewidth',1.0,'fontsize',11,'fontname','Arial');
[caz,cel] = view()
%view([-37.5 20])
%%
RESIRE_Dose1p6e5=load('Resire_aSi2nm_FosCetr_Dose1p6e5_cov35pm_BM3D_S0p25N200.mat')
rec_Dose1p6e5=RESIRE_Dose1p6e5.OBJ.reconstruction;
mask=load('Mask_aSi_cov35pm.mat');
mask=mask.mask_BB;
rec_Dose1p6e5=mask.*rec_Dose1p6e5;

%figure('Renderer', 'painters', 'Position', [100 100 300 300]); 
figure
clf; 
set(gcf,'color','w');
set(gcf,'outerposition',[100 100 300 300])
isosurface(rec_Dose1p6e5,0.35)
axis image;
xlim([64-50,65+50]); ylim([64-50,65+50]); zlim([64-50,65+55]);
xticks(20:20:100); yticks(20:20:100); zticks(20:20:100);
colormap hot
axis tight
axis equal
ax = gca;
ax.LineWidth = 1.0;
set(gca,'linewidth',1.0,'fontsize',11,'fontname','Arial');
%view([-50 15])
%%
Proj_Dose8p0e3_Noise=load('aSi_2nmFosCetr_Dose8p0e3_cov35pm_ForBM3D.mat');
Proj_Dose8p0e3_Noise=Proj_Dose8p0e3_Noise.Trunc_projections_forBM3D;
Proj_Dose8p0e3_BM3D=RESIRE_Dose8p0e3.OBJ.InputProjections;
Proj_Dose8p0e3_Cal=RESIRE_Dose8p0e3.OBJ.Calprojection;

figure('Renderer', 'painters', 'Position', [100 100 600 220])
img(Proj_Dose8p0e3_Noise,[],Proj_Dose8p0e3_BM3D,[],Proj_Dose8p0e3_Cal,[],'colormap','gray')

%%
Proj_Dose1p6e4_Noise=load('aSi_2nmFosCetr_Dose1p6e4_cov35pm_ForBM3D.mat');
Proj_Dose1p6e4_Noise=Proj_Dose1p6e4_Noise.Trunc_projections_forBM3D;
Proj_Dose1p6e4_BM3D=RESIRE_Dose1p6e4.OBJ.InputProjections;
Proj_Dose1p6e4_Cal=RESIRE_Dose1p6e4.OBJ.Calprojection;

figure('Renderer', 'painters', 'Position', [100 100 600 220])
img_plot(Proj_Dose1p6e4_Noise,[],Proj_Dose1p6e4_BM3D,[],Proj_Dose1p6e4_Cal,[],'colormap','gray')
%print('-r500','-dpdf','si_multislice_defocus_0_8e4_ct_iso_0_3.pdf');
print('-dbmp', '-r600', 'Fig_1a-c.bmp');
%%
Proj_Dose1p6e5_Noise=load('aSi_2nmFosCetr_Dose1p6e5_cov35pm_ForBM3D.mat');
Proj_Dose1p6e5_Noise=Proj_Dose1p6e5_Noise.Trunc_projections_forBM3D;
Proj_Dose1p6e5_BM3D=RESIRE_Dose1p6e5.OBJ.InputProjections;
Proj_Dose1p6e5_Cal=RESIRE_Dose1p6e5.OBJ.Calprojection;

figure('Renderer', 'painters', 'Position', [100 100 600 220])
img(Proj_Dose1p6e5_Noise,[],Proj_Dose1p6e5_BM3D,[],Proj_Dose1p6e5_Cal,[],'colormap','gray')




%%
