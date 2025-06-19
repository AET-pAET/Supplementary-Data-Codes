% this code calculate shift from the center of mass method
close all; clear all;
%%
addpath('../calculate_shiftXY');
addpath('./src');
%%
imgStackDeblur=load('Proj_imgStackDeblur_VST');
imgStackDeblur=imgStackDeblur.imgStackDeblur;
imgStackDeblur=imgStackDeblur(:,:,1:55);
%%
COMrefX = round((size(imgStackDeblur,1)+1)/2);
COMrefY = round((size(imgStackDeblur,2)+1)/2);
CoM_Loct = zeros(2,55);

imgStackDeblur_shift = zeros(size(imgStackDeblur));
for i = 1:size(imgStackDeblur,3)
    [COMX, COMY] = My_COM(squeeze(imgStackDeblur(:,:,i)));

    xshift = COMrefX - COMX;
    yshift = COMrefY - COMY;
    CoM_Loct(1:2,i)=[xshift yshift];

    imgStackDeblur_shift(:,:,i) = real(My_FourierShift(imgStackDeblur(:,:,i),xshift,yshift));
end
%%
%save('CoPdPt_test_align_COM_data.mat','CoM_Loct');
%%
center_corr=[0.25 -0.1];% center correction from GT.

figure
plot(-CoM_Loct(1,:)*0.333+center_corr(1))
hold on
plot(-CoM_Loct(2,:)*0.333+center_corr(2))
%% compare shift from COM and Model
load('CoPdPt_test_align_data.mat', 'shiftX', 'shiftY'); %Alignment by Center of Mass
figure(11); clf; set(gcf,'position',[250,250,300,250]);
hold on;
plot(shiftX*0.333); plot(shiftY*0.333); box on;

plot(-CoM_Loct(1,:)*0.333+center_corr(1))
hold on
plot(-CoM_Loct(2,:)*0.333+center_corr(2))

legend({'ShiftX','ShiftY'})

xlim([0,55]); ylim([-1.2,2.2]); yticks(-1:1);
legend({'ShiftX','ShiftY'})
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(gca, 'Layer', 'top')
legend boxoff   
xlabel(['Image numbers'])
ylabel('Shift (Å)')