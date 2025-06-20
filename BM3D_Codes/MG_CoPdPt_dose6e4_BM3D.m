clear; clc;
%%
addpath('./BM3D')
%%
Proj=load('CoPdPt_imgStack.mat');
Proj_Dose5p6e4=Proj.imgStack;
% figure
% img(Proj_Dose5p6e4(1:end,:,1:end),'colormap','gray')
%%
CropDataAll_reduDarkC=double(Proj_Dose5p6e4(:,:,:));
[Alpha_arr,Sigma_arr] = Main_getAlphaSigma_parameters(CropDataAll_reduDarkC);
Trunc_projections_forBM3D = CropDataAll_reduDarkC;

Trunc_projections_forBM3D =Trunc_projections_forBM3D;
%%
alpha_mean = mean(Alpha_arr);
sigma_mean = mean(Sigma_arr);
Dset_AfterBM3D_bg = BM3D_Main(Trunc_projections_forBM3D(:,:,:),1,alpha_mean*ones(size(Alpha_arr))*0.55,0.55*sigma_mean*ones(size(Sigma_arr)));
%%

% save('CoPdPt_Proj_Dose5p6e4_BM3D_NoCOM.mat','Dset_AfterBM3D_bg')


%% Here we show an example of BM3D using our optimized parameters. Other datasets can do in the same way.






