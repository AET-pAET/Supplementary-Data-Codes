
% Data2_tomo1 = Projections_NiPt_1118_tomo1_AfterImgReg_Trunc_400;
function [Alpha_arr,Sigma_arr,varargout] = Main_getAlphaSigma_parameters(Data2,varargin)
MaskStack = zeros(size(Data2,1),size(Data2,2),size(Data2,3));

if nargin > 2
    Alpha_name = varargin{1};
    Sigma_name = varargin{2};
    savefile = 1;
else
    savefile = 0;
end
Alpha_arr = zeros(1,size(Data2,3));
Sigma_arr = zeros(1,size(Data2,3));
% figure;
for pj = 1:size(Data2,3)
    disp(pj);
    Data = squeeze(Data2(:,:,pj,:));
    % create mask
    A1=edge_preserve_smoothing_YY(squeeze(Data(:,:,1)),3,1);
    % inital Otsu threshold 
    im=double(A1/max(max(A1))*255);  % convert to 0:255 range to fit into otsu method
    [ot,~]=otsu_thresh(im);
    ot=ot*max(max(A1))/255*0.9; % use loose otsu threshold to choose far away region
    after_conv = (A1>ot) * 1;
    after_conv  = bwareaopen(after_conv, 3000);   %%% there are also circularity, long to short radius ratio, etc.
    se = strel('disk',20);
    after_conv = imdilate(after_conv,se);
    after_conv = bg_conv_gauss_corr(after_conv, 3, 0.1); %loose mask
%     after_conv = zeros(size(Data,1),size(Data,2));
    currMask = ~after_conv;
%     currMask = ones(size(Data,1),size(Data,2)) - after_conv;
    MaskStack(:,:,pj) = currMask;
    
%     imagesc(A1);axis image;title(sprintf('initial %d',pj));pause
%     imagesc(A1.*currMask);axis image;title(sprintf('initial %d',pj));pause

    AlpBinSize = 2;
    AlpThFact  = 1;
    SigBinSize = 2;
    [Alpha, Sigma, curr_mean] = My_FePt_Parameter_Estimation_BGmask_varAnalYY(Data,AlpBinSize,AlpThFact,SigBinSize,currMask);
    fprintf('Alpha is %d \n',Alpha)
    fprintf('Sigma is %d \n',Sigma)
    Alpha_arr(pj) = Alpha;
    Sigma_arr(pj) = Sigma;
     
end
%%
if nargout >2
    varargout{1} = MaskStack;
end
%%
if savefile == 1
save(Alpha_name, 'Alpha_arr', '-mat');
save(Sigma_name, 'Sigma_arr', '-mat');
end
end