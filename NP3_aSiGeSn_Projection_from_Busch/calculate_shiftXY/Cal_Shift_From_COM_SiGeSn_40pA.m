close all; clear all;
%%
load('SiGeSn_40pA_imgStackDeblur.mat','imgStackDeblur');
%%
COMrefX = round((size(imgStackDeblur,1)+1)/2);
COMrefY = round((size(imgStackDeblur,2)+1)/2);
CoM_Loct = zeros(2,60);

imgStackDeblur_shift = zeros(size(imgStackDeblur));
for i = 1:size(imgStackDeblur,3)
    [COMX, COMY] = My_COM(squeeze(imgStackDeblur(:,:,i)));

    xshift = COMrefX - COMX;
    yshift = COMrefY - COMY;
    CoM_Loct(1:2,i)=[xshift yshift];

    imgStackDeblur_shift(:,:,i) = real(My_FourierShift(imgStackDeblur(:,:,i),xshift,yshift));
end
%%
save('SiGeSn_40pA_align_2_data.mat','CoM_Loct');







