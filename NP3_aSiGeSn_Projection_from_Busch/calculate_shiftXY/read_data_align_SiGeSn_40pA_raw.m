fclose('all');
num_img = 60;
name = 'Test_Particle';
img_size = 512;

filepath = './SiGeSn-1.6E4e_per_A2';
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
%%
% imgStack=single(imgStack);
% save('SiGeSn_40pA_imgStack.mat','imgStack');
%%
COMrefX = round((size(imgStack,1)+1)/2);
COMrefY = round((size(imgStack,2)+1)/2);
CoM_Loct = zeros(2,60);

imgStack_shift = zeros(size(imgStack));
for i = 1:size(imgStack,3)
    [COMX, COMY] = My_COM(squeeze(imgStack(:,:,i)));

    xshift = COMrefX - COMX;
    yshift = COMrefY - COMY;
    CoM_Loct(1:2,i)=[xshift yshift];

    imgStack_shift(:,:,i) = real(My_FourierShift(imgStack(:,:,i),xshift,yshift));
end

% save('SiGeSn_40pA_align_1.mat','imgStack_shift', 'CoM_Loct');
%%
imgStackDeblur = zeros(size(imgStack));
for k = 1:num_img
    k
    baseFileName = ['rotation_',num2str(k-1)];
    imageArray = imgStack(:,:,k);%imgStack(:,:,k);
    %imageArray = Up_imgStack(:,:,k);

    % load image to restore
    y=imageArray;

    % choose image peak value
    peak = 255;
    y(y<=0)=0;
    y=y/max(y(:))*peak;
    % blur point-spread function (PSF) and blurring
    v = fspecial('gaussian', 25, sqrt(9)); %25 original
    g = imfilter(y, v(end:-1:1,end:-1:1), 'circular');

    % seed for noise realization
    ii=1; randn('seed',ii);  rand('seed',ii);
    z=poissrnd(g);

    % Deblur via iterative (VST)+(colored noise removal)
    % yhat_RWI is the final estimate
    % yhat_RI is the estimate after the first stage (standard regularized inverse)
    [yhat_RWI,yhat_RI] = iterVSTpoissonDeb(z,v);

    %  show also intermediate RI result, or only RWI
    showRI=1;

    % compute PSNR of the estimate
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
    % filtered_image = BPF(yhat_RWI,1,128,3);
    filtered_image = yhat_RWI;%BPF(yhat_RWI,1,128,3);
    filtered_image = double(filtered_image);
    filtered_image = rescale(filtered_image,0,1);
    filtered_image = imgaussfilt(imageArray,1);
    imgStackDeblur(:,:,k)  = rescale(filtered_image,0,1);
%     saveFileName = ['output/deblur/',baseFileName,'.tiff'];
%     imwrite(imgStackDeblur(:,:,k),saveFileName,'tiff'); %tiff output will be integer variable maintains double
    imshow(imgStackDeblur(:,:,k));
    drawnow;
end
%%
% imgStackDeblur=single(imgStackDeblur);
% save('SiGeSn_40pA_imgStackDeblur.mat','imgStackDeblur');
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
%save('SiGeSn_40pA_align_2.mat','CoM_Loct');

