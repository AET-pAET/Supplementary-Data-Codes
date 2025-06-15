%%%% OpTomo Control file - Robert Busch University of Illinois 2024-12-05
%% This Setup requires astra-toolbox and opSpot
%% Astra-toolbox install at: https://github.com/astra-toolbox/astra-toolbox
%% modified from https://github.com/astra-toolbox/astra-toolbox/blob/master/matlab/tools/opTomo.m
%% See PAPER: https://doi.org/10.1007/s11075-015-0016-4
%% opSpot http://www.cs.ubc.ca/labs/scl/spot/
%% Setup
filepath = 'C:\Path\To\file';
name = 'Test_Particle';
num_img = 55;
img_size = 300;
up_rate = 1;
angles = deg2rad(angles);
maxit = 200;

%%
if ~isdir(filepath)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', filepath);
  uiwait(warndlg(errorMessage));
  return;
end

%% Load File
fclose('all');

for k = 1:num_img
    k
  baseFileName = ['rotation_',num2str((k-1)),'.img'];
  fullFileName = fullfile(filepath, baseFileName);
  fileID = fopen(fullFileName);  
  imageArray = fread(fileID,[img_size img_size],'float');

  imshow(imageArray);  % Display image.
  imageArray(imageArray<0.0)=0;
  drawnow; % Force display to update immediately.
  
  imgStack(:,:,k) = imageArray;
  
  fclose(fileID);
end

max(imgStack(:))
save('imgStack.mat','imgStack');
%clear all
%%
clear imgStackDeblur
mkdir output
mkdir output/deblur


%% Image Denois IVST-3BMD
for k = 1:num_img
  k
  baseFileName = ['rotation_',num2str(k-1)];
  imageArray = imgStack(:,:,k);%imgStack(:,:,k);
  %imageArray = Up_imgStack(:,:,k);
  
%% load image to restore 
    y=imageArray;

%% choose image peak value
    peak = 255;
    y(y<=0)=0;
    y=y/max(y(:))*peak;

%% blur point-spread function (PSF) and blurring
    v = fspecial('gaussian', 25, sqrt(9)); %25 original
    g = imfilter(y, v(end:-1:1,end:-1:1), 'circular'); 

%% seed for noise realization
    ii=1; randn('seed',ii);  rand('seed',ii);
    z=poissrnd(g);

%% Deblur via iterative (VST)+(colored noise removal)
% yhat_RWI is the final estimate
% yhat_RI is the estimate after the first stage (standard regularized inverse)
    [yhat_RWI,yhat_RI] = iterVSTpoissonDeb(z,v);

%%  show also intermediate RI result, or only RWI
    showRI=1;

%% compute PSNR of the estimate
    PSNR_z=10*log10(peak^2/mean((y(:)-z(:)).^2));
    PSNR_yhat_RI=10*log10(peak^2/mean((y(:)-yhat_RI(:)).^2));
    PSNR_yhat_RWI=10*log10(peak^2/mean((y(:)-yhat_RWI(:)).^2));
    resString=['PSNR(z): '     num2str(PSNR_z,'%4.2f'),'dB'];
    if showRI
        resString=[resString, '   PSNR(yhat_RI): '  num2str(PSNR_yhat_RI,'%4.2f'),'dB'];
    end
    resString=[resString, '   PSNR(yhat_RWI): ' num2str(PSNR_yhat_RWI,'%4.2f'),'dB'];
    disp(resString);
    clf
    filtered_image = BPF(yhat_RWI,1,128,3);
    filtered_image = double(filtered_image);
    filtered_image = rescale(filtered_image,0,1);
    filtered_image = imgaussfilt(imageArray,1);
    imgStackDeblur(:,:,k)  = rescale(filtered_image,0,1);
    saveFileName = ['output/deblur/',baseFileName,'.tiff'];
    imwrite(imgStackDeblur(:,:,k),saveFileName,'tiff'); %tiff output will be integer variable maintains double
    imshow(imgStackDeblur(:,:,k));
    drawnow;
    
end
save('imgStackDeblur.mat','imgStackDeblur');

%% Interpolation of Denoised images
clear X Y XX YY Up_imgStack Upsample
%Upsample = imgStack;
Upsample = imgStackDeblur;
X = 1:size(Upsample,1);
Y = 1:size(Upsample,2);
x_up = up_rate;               %upsample values
y_up = up_rate;            %upsample values
imgsize = size(Upsample,3);
XX = linspace(1,size(Upsample,1),(size(Upsample,1)*x_up));
YY = linspace(1,size(Upsample,2),(size(Upsample,2)*y_up));
Up_imgStack = zeros(size(Upsample,1)*x_up, size(Upsample,2)*y_up, imgsize);
I1 = zeros(size(Upsample,1),(size(Upsample,2)*y_up));
I2 = zeros(size(Upsample,1)*x_up,(size(Upsample,2)*y_up));
for i = 1:imgsize
    I = Upsample(:,:,i);
    for j = 1:size(Upsample,1)
        I1(j,:) = spline(X, I(j,:), XX);
    end
    for k = 1:size(Upsample,1)*y_up
        I2(:,k) = spline(Y, I1(:,k), YY);
    end
    Up_imgStack(:,:,i) = I2;
end
clear X Y XX YY  Upsample y_up x_up I I1 I2 i j k imgsize
save('Up_imgStack.mat')

%% Setup geometry for SIRT
%p = permute(imgStackDeblur, [2 3 1]);
p = permute(Up_imgStack, [2 3 1]); % 1 3 2, 2 3 1

p = rescale(p,0,1);
p(p <= 0) = 0;
p = rescale(p,0,1);


l = numel(angles);
vectors = zeros(l,12);

%y axis tilt, Currently these are written for single axis rotation

for i = 1:l
    a = angles(i);
    vectors(i,1:3) = [sin(a), 0, -cos(a)];
    vectors(i,4:6) = [0, 0, 0];
    vectors(i,7:9) = [cos(a), 0, sin(a)];
    vectors(i,10:12) = [0, 1, 0];
end

%x axis tilt
%for i = 1:l
%    a = angles(i);
%    vectors(i,1:3) = [0, -sin(a), -cos(a)];
%    vectors(i,4:6) = [0, 0, 0];
%    vectors(i,7:9) = [0, -cos(a),  sin(a)];
%    vectors(i,10:12) = [1, 0, 0];
%end



d = size(p, 1);
proj_geom = astra_create_proj_geom('parallel3d_vec', d, d, vectors);

s = d;
vol_geom = astra_create_vol_geom([s s s]);

vol_id = astra_mex_data3d('create', '-vol', vol_geom, 0);
proj_id = astra_mex_data3d('create', '-proj3d', proj_geom, p);


%% modified from https://github.com/astra-toolbox/astra-toolbox/blob/master/matlab/tools/opTomo.m
%% See PAPER: https://doi.org/10.1007/s11075-015-0016-4
%% opSpot http://www.cs.ubc.ca/labs/scl/spot/
% These operators (C in particular) can be memory intensive, once
% established all memory for the calculations is reserved
clear imgStack imgStackDeblur 
W = opTomo('cuda', proj_geom, vol_geom);

c = 1./sum(W,1);
c(c==Inf) = 0;
C = opDiag(c);
clear c


r = 1./sum(W,2);
r(r==Inf) = 0;
R = opDiag(r);
clear r

%%
p = p(:);
v = zeros(d,d,d,'double');
v = v(:);
for i = 1:maxit
    v = v + C*W'*R*(p - W*v); 
    v(v<=0) = 0; %Negativity constraint if wanted. This breaks garunteed  
                 %convergence but is generally regarded as giving more 
                 %accurate reconstructions
    i
end

SIRT_int = reshape(v, W.vol_size);
%save(['output/',name,'.mat'],'SIRT');
%%
[sinogram_id, sinogram] = astra_create_sino3d_cuda(SIRT_int, proj_geom, vol_geom);

%%
SIRT = imresize3(SIRT_int,[img_size img_size img_size],'cubic');
%%
%save(['output/',name,'.mat'],'SIRT');
%%
%save(['output/',name,'_sinogram.mat'],'sinogram');
%%
%clearvars -except angles SIRT_int imgStackDeblur sinogram imgStackDeblur SIRT
%% From Here, We run the polynomial tracing from Yang et. al.
