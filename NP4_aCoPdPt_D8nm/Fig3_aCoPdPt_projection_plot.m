close all; clear all;
%%
addpath('4.1_projections');
addpath('4.2_reconstructions');
%%
RESIRE_Dose0=load('Resire_CoPdPt_NoRadmDefocus_cov35pm_Dose1p7e4_BM3DOPT_S1p0N1000.mat')
Proj_Dose0_Noise=load('Proj_CoPdPt_NoRadmDefocus_cov35pm_Dose1p7e4_BeforeBM3D.mat');
Proj_Dose0_Noise=Proj_Dose0_Noise.Proj_Dose1p7e4_PoissonNoise;
%%
Proj_Dose1p7e4_BM3D=RESIRE_Dose0.OBJ.InputProjections;
Proj_Dose1p7e4_Cal=RESIRE_Dose0.OBJ.Calprojection;

figure('Renderer', 'painters', 'Position', [100 100 600 220])
img_plot(Proj_Dose0_Noise,[],Proj_Dose1p7e4_BM3D,[],Proj_Dose1p7e4_Cal,[],'colormap','gray')

print('-dbmp', '-r1200', 'Fig_3a-c.bmp');
print('-djpeg', '-r1200', 'Fig_3a-c.jpg');




return
