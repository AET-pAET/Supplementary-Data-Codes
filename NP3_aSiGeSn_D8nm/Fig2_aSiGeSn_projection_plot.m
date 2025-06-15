close all; clear all;
%%
addpath('3.1_projections');
addpath('3.2_reconstructions');
%%
RESIRE_Dose1p6e4=load('Resire_aSiGeSn_FosCetr_Dose1p6e4_cov35pm_BM3D_S2p0N2000.mat')
rec_Dose1p6e4=RESIRE_Dose1p6e4.OBJ.reconstruction;
%%
RESIRE_Dose1p6e5=load('Resire_aSiGeSn_FosCetr_Dose1p6e5_cov35pm_BM3D_S2p0N2000.mat')
rec_Dose1p6e5=RESIRE_Dose1p6e5.OBJ.reconstruction;
%%
Proj_Dose1p6e4_Noise=load('aSiGeSn_FosCetr_Dose1p6e4_cov35pm_ForBM3D.mat');
Proj_Dose1p6e4_Noise=Proj_Dose1p6e4_Noise.Trunc_projections_forBM3D;
Proj_Dose1p6e4_BM3D=RESIRE_Dose1p6e4.OBJ.InputProjections;
Proj_Dose1p6e4_Cal=RESIRE_Dose1p6e4.OBJ.Calprojection;

figure('Renderer', 'painters', 'Position', [100 100 600 220])
img_plot(Proj_Dose1p6e4_Noise,[],Proj_Dose1p6e4_BM3D,[],Proj_Dose1p6e4_Cal,[],'colormap','gray')

%print('-dbmp', '-r600', 'Fig_2a-c.bmp');
print('-djpeg', '-r1200', 'Fig_2a-c.jpg');
%%
Proj_Dose1p6e5_Noise=load('aSiGeSn_FosCetr_Dose1p6e5_cov35pm_ForBM3D.mat');
Proj_Dose1p6e5_Noise=Proj_Dose1p6e5_Noise.Trunc_projections_forBM3D;
Proj_Dose1p6e5_BM3D=RESIRE_Dose1p6e5.OBJ.InputProjections;
Proj_Dose1p6e5_Cal=RESIRE_Dose1p6e5.OBJ.Calprojection;

figure('Renderer', 'painters', 'Position', [100 100 600 220])
img(Proj_Dose1p6e5_Noise,[],Proj_Dose1p6e5_BM3D,[],Proj_Dose1p6e5_Cal,[],'colormap','gray')



