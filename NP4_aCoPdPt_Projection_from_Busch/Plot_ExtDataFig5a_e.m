close all; clear all
%%
RESIRE_Data=load('3_Resire_Dose5p6e4_BM3D_After_Alignment.mat')
RESIRE_Rec=RESIRE_Data.OBJ.reconstruction;

SIRT_Rec = importdata('MG_5p6e4_SIRT_Reconstruction_Provided_by_Busch.mat');
SIRT_Rec = flip(rot90(SIRT_Rec,1),1);
%%
img(sum(SIRT_Rec(:,:,151:151),3),[],sum(flip(rot90(RESIRE_Rec(:,:,151:151))),3),[],'colormap','gray');
%print('-r600','-djpeg','ExtData_Fig5de.jpg');
%%
load('CoPdPt_test_align.mat', 'shiftX', 'shiftY'); %Alignment by Center of Mass
figure(11); clf; set(gcf,'position',[250,250,300,250]);
hold on;
plot(shiftX*0.347); plot(shiftY*0.347); box on;
legend({'ShiftX','ShiftY'})

xlim([0,55]); ylim([-1.2,2.2]); yticks(-1:1);
legend({'ShiftX','ShiftY'})
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(gca, 'Layer', 'top')
legend boxoff   
xlabel(['Image numbers'])
ylabel('Shift (Å)')
%print('-r600','-djpeg','ExtData_Fig5a.jpg');
%%
load('CoPdPt_Misaligned_Denoised_Projection_using_Busch_Code.mat')
plotCommonLine(imgStackDeblur,2)
set(gcf,'position',[250,250,300,250]);
xlim([0,300]);
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(gca, 'Layer', 'top')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
xlabel(['Y (pixel)'])
ylabel('Summed Intensity')
xticks(0:100:500)
%print('-r600','-djpeg','ExtData_Fig5b.jpg');
%%
Projection_After_Alignment=RESIRE_Data.OBJ.InputProjections/500;

for i=1:55
Projection_After_Alignment(:,:,i)=Projection_After_Alignment(:,:,i)-min(min(Projection_After_Alignment(:,:,i)));
end

plotCommonLine(Projection_After_Alignment,1)

set(gcf,'position',[250,250,300,250]);
xlim([0,300]); ylim([0,90]);
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(gca, 'Layer', 'top')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
xlabel(['Y (pixel)'])
ylabel('Summed Intensity')
xticks(0:100:300)

%print('-r600','-djpeg','ExtData_Fig5c.jpg');

