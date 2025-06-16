close all; clear all;
%%
addpath('5.1_projections');
addpath('5.2_reconstructions');
%%
proj_single=load('phases_single_slice_1.6e3_SiO2_cov35pm.mat');
proj_single=proj_single.projections_cov35pm;

proj_ml=load('phases_multislice_1.6e3_SiO2_cov35pm.mat');
proj_ml=proj_ml.projections_cov35pm;


figure('Renderer', 'painters', 'Position', [100 100 600 200])
img(proj_ml(1:end,1:end,:),[],proj_single(1:end,1:end,:),[],'colormap','gray')
%%
figure('Renderer', 'painters', 'Position', [100 100 400 150])
img_plot(proj_ml(1:end,1:end,:),[],proj_single(1:end,1:end,:),[],'colormap','gray')

print('-dbmp', '-r1200', 'Fig_4a-b.bmp');
print('-djpeg', '-r2400', 'Fig_4a-b.jpg');