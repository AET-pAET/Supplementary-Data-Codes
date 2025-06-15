close all; clear all
%%
load('SiGeSn_40pA_reconres_busch_denoised_aligned.mat', 'rec3');% Reconstruction after alignment using RESIRE
load('SiGeSn_40pA_10us.mat');% Reconstruction before alignment using Busch's projection and code
%% compare reconstruction
%SIRT_rot90 = flip(rot90(SIRT,1),1);
SIRT_int_rot90 = flip(rot90(SIRT_int,1),1);
img(sum(SIRT_int_rot90(:,:,255:259),3),[],sum(rec3(:,:,255:259),3),[],'colormap','gray');
%
%print('-r600','-djpeg','ExtData_Fig4de.jpg');
%% plot CM
load('SiGeSn_40pA_align_2.mat', 'CoM_Loct')
figure(11); set(gcf,'position',[250,250,300,250]);
clf; hold on;
plot(-CoM_Loct'*0.195); box on;
xlim([0,61]); ylim([-1.5,1.5]); yticks(-1:1);
legend({'ShiftX','ShiftY'})
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(gca, 'Layer', 'top')
legend boxoff   
xlabel(['Image numbers'])
ylabel('Shift (Å)')
%print('-r600','-djpeg','ExtData_Fig4a.jpg');
%% Busch_commonline
load('SiGeSn_40pA_imgStackDeblur.mat')
plotCommonLine(imgStackDeblur,2)
set(gcf,'position',[250,250,300,250]);
xlim([0,513]);
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(gca, 'Layer', 'top')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
xlabel(['Y (pixel)'])
ylabel('Summed Intensity')
xticks(0:100:500)
%print('-r600','-djpeg','ExtData_Fig4b.jpg');
%% commonline after alignment
load('SiGeSn_40pA_bm3d_mask_shift2.mat')
plotCommonLine(Dset_ys_mask_shift2,2)
set(gcf,'position',[250,250,300,250]);
xlim([0,513]); ylim([0,14000]);
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(gca, 'Layer', 'top')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
xlabel(['Y (pixel)'])
ylabel('Summed Intensity')
xticks(0:100:500)
%print('-r600','-djpeg','ExtData_Fig4c.jpg');

