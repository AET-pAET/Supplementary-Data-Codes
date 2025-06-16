close all; clear all;
%%
addpath('2.1_projections');
addpath('2.2_reconstructions');
%%
RESIRE_Dose8p0e3=load('Resire_aSi9nm_FosCetr_Dose8p0e3_cov35pm_BM3D_S1p0N200.mat')
rec_Dose8p0e3=RESIRE_Dose8p0e3.OBJ.reconstruction;
%%
RESIRE_Dose1p6e4=load('Resire_aSi9nm_FosCetr_Dose1p6e4_cov35pm_BM3D_S1p0N200.mat')
rec_Dose1p6e4=RESIRE_Dose1p6e4.OBJ.reconstruction;
%%
RESIRE_Dose1p6e5=load('Resire_aSi9nm_FosCetr_Dose1p6e5_cov35pm_BM3D_S1p0N200.mat')
rec_Dose1p6e5=RESIRE_Dose1p6e5.OBJ.reconstruction;
%%
Proj_Dose8p0e3_Noise=load('aSi_9nmFosCetr_Dose8p0e3_cov35pm_ForBM3D.mat');
Proj_Dose8p0e3_Noise=Proj_Dose8p0e3_Noise.Trunc_projections_forBM3D;
Proj_Dose8p0e3_BM3D=RESIRE_Dose8p0e3.OBJ.InputProjections;
Proj_Dose8p0e3_Cal=RESIRE_Dose8p0e3.OBJ.Calprojection;

figure('Renderer', 'painters', 'Position', [100 100 600 220])
img_plot(Proj_Dose8p0e3_Noise,[],Proj_Dose8p0e3_BM3D,[],Proj_Dose8p0e3_Cal,[],'colormap','gray')

%%
Proj_Dose1p6e4_Noise=load('aSi_9nmFosCetr_Dose1p6e4_cov35pm_ForBM3D.mat');
Proj_Dose1p6e4_Noise=Proj_Dose1p6e4_Noise.Trunc_projections_forBM3D;
Proj_Dose1p6e4_BM3D=RESIRE_Dose1p6e4.OBJ.InputProjections;
Proj_Dose1p6e4_Cal=RESIRE_Dose1p6e4.OBJ.Calprojection;

figure('Renderer', 'painters', 'Position', [100 100 600 220])
img(Proj_Dose1p6e4_Noise,[],Proj_Dose1p6e4_BM3D,[],Proj_Dose1p6e4_Cal,[],'colormap','gray')
print('-dbmp', '-r600', 'Fig_1e-g.bmp');
%%
Proj_Dose1p6e5_Noise=load('aSi_9nmFosCetr_Dose1p6e5_cov35pm_ForBM3D.mat');
Proj_Dose1p6e5_Noise=Proj_Dose1p6e5_Noise.Trunc_projections_forBM3D;
Proj_Dose1p6e5_BM3D=RESIRE_Dose1p6e5.OBJ.InputProjections;
Proj_Dose1p6e5_Cal=RESIRE_Dose1p6e5.OBJ.Calprojection;

figure('Renderer', 'painters', 'Position', [100 100 500 220])
img(Proj_Dose1p6e5_Noise,[],Proj_Dose1p6e5_BM3D,[],Proj_Dose1p6e5_Cal,[],'colormap','gray')