close all; clear all;
%%
addpath('1.11_projections');
addpath('1.12_reconstructions');
%%
Proj_raw_Busch=load("Proj_aSi_40pA_raw_Busch.mat");
Proj_raw_Busch=Proj_raw_Busch.imgStack;

figure
img(Proj_raw_Busch(1:end,:,:),[],'colormap','gray')

for i=1:60
Proj_raw(:,:,i)=Proj_raw_Busch(:,:,i)';
end
% We corped the nanopartilce region for BM3D and reconstruction
Proj_Dose40pA_raw=Proj_raw(512/2-64+1:512/2+64,512/2-64+1:512/2+64,:); 
Proj_Dose40pA_raw=Proj_Dose40pA_raw-0;
Proj_Dose40pA_raw(Proj_Dose40pA_raw<0)=0;

figure
img(Proj_Dose40pA_raw(1:end,:,:),[],'colormap','gray')
%%
RESIRE_Dose40pA=load("Resire_Si_40pA_BM3D.mat");
Proj_Dose40pA_BM3D=RESIRE_Dose40pA.OBJ.InputProjections;
Proj_Dose40pA_Cal=RESIRE_Dose40pA.OBJ.Calprojection;

figure('Renderer', 'painters', 'Position', [100 100 600 220])
img_plot(Proj_Dose40pA_raw,[],Proj_Dose40pA_BM3D,[],Proj_Dose40pA_Cal,[],'colormap','gray')