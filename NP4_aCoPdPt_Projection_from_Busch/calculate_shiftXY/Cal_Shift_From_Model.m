clear all; close all;
%%
addpath('./src')
%%
% model_OriOri = importdata('model_refined_HEA_181027_t1_yk55_yk_0829.mat');
% model_OriOri = double(model_OriOri);
% temp_atomtype = importdata('LocalC_HEA_181027_yk55_yk_gd_localtype.mat');
% 
% euler_angles = zeros(70,3);
% euler_angles(:,3) = 0:2.6:180;
% xdata.angles = euler_angles;
% xdata.model = model_OriOri;
% xdata.atoms = temp_atomtype;
% xdata.Z_arr = [27,46,78];
% xdata.Res = 0.333;
% xdata.volSize = 300;
% xdata.CropHalfWidth = 3;
% 
% para = zeros(2,numel(temp_atomtype));
% para(1,temp_atomtype==1) = 0.5;
% para(1,temp_atomtype==2) = 1.0;
% para(1,temp_atomtype==3) = 1.5;
% para(2,temp_atomtype==1) = 3;
% para(2,temp_atomtype==2) = 3;
% para(2,temp_atomtype==3) = 3;
% 
% Projs = calcLinearPjs_fit_para(para,xdata);
% %
% Projs=single(Projs);
% save('Projs_from_GT.mat',"Projs");
%%
% Obtain the Projs by running the above code.
Projs=load('Projs_from_GT.mat');
Projs=Projs.Projs;
% Obtain the imgStack by running the Busch's code in folder "SIRT_Code_Used_by_Busch."
imgStack=load('imgStack.mat'); 
imgStack=imgStack.imgStack;
% Obtain the imgStackDeblur by running the Busch's code in folder "SIRT_Code_Used_by_Busch."
imgStackDeblur=load('Proj_imgStackDeblur_VST.mat');
imgStackDeblur=imgStackDeblur.imgStackDeblur;
%%
figure(11); clf;
for i = 1:55
    imagesc(Projs(:,:,i)); axis image;
    title(sprintf('simulated pj %d',i)); pause();
    
    imagesc(imgStack(:,:,i)); axis image;
    title(sprintf('Busch''s pj %d',i)); pause();
end
%% obtain shift by cross-correlation
S = imgStack;
cen = round((size(S,1)+1)/2);
metric_arr=ones(1,55);
shiftX = zeros(1,55);
shiftY = zeros(1,55);
for i = 1:55
    disp(num2str(i))
    [metric, suggested_center_x, suggested_center_y] = alignByNormXCorrSubpixel(S(:,:,i),Projs(:,:,i),0.05);
    metric_arr(i)=metric;
    shiftX(i) = suggested_center_x-cen;
    shiftY(i) = suggested_center_y-cen;
    
    S(:,:,i) = real(My_RealShift(S(:,:,i),-shiftX(i), -shiftY(i)));
end
%% obtain shift by cross-correlation
S_2 = imgStackDeblur;
cen = round((size(S_2,1)+1)/2);
metric_arr=ones(1,55);
shiftX = zeros(1,55);
shiftY = zeros(1,55);
for i = 1:55
    disp(num2str(i))
    [metric, suggested_center_x, suggested_center_y] = alignByNormXCorrSubpixel(S_2(:,:,i),Projs(:,:,i),0.05);
    metric_arr(i)=metric;
    shiftX(i) = suggested_center_x-cen;
    shiftY(i) = suggested_center_y-cen;
    
    S_2(:,:,i) = real(My_RealShift(S_2(:,:,i),-shiftX(i), -shiftY(i)));
end
%% save the data

%save('CoPdPt_test_align_data.mat','shiftX', 'shiftY')
%%
figure
plot(shiftX*0.333)
hold on
plot(shiftY*0.333)



