function oim=edge_preserve_smoothing_YY(im,ws,flat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% oim=edge_preserve_smoothing(im,ws,flat)
%   where im is the image (grey)
%       ws is the window size (2ws+1)x(2ws+1) (default, ws=2)
%       flat=1 indicates use of flat smoothing instead of Gaussian
%           Edge preserving smoothing is performed on the input image
%       oim is the smoothed image

if nargin<2
    ws=2;
end

if nargin<3
    flat=0;
end

[height,width]=size(im);
ww=1+2*ws; %window width
wp=ww^2; %number of window pixels

%find the variance of the image in every ww by ww window
varnim=conv2(im.^2,ones(ww)/wp,'same')-conv2(im,ones(ww)/wp,'same').^2;

if flat~=1
    gfilter=ch4_gaussian_filter(ws,0.1); %find a Gaussian filter ws by ws with max edge cutoff of 0.1
end

oim=im; %in the output image there will be a border of unproccessed pixels all around

height_cut = height-2*ww+2;
width_cut = width-2*ww+2;
varover = zeros(height_cut,width_cut,wp);

height_diff = height - height_cut;
width_diff = width - width_cut;
x_add = height_diff/2 - ws;
y_add = width_diff/2 - ws;
for i=1:ww
    for j=1:ww
        varover(:,:,i+(j-1)*ww) = varnim(i+x_add:i+x_add+height_cut-1,j+y_add:j+y_add+width_cut-1);
    end
end

segment_ind = zeros(height_cut,width_cut,2);
[~, vi] = min(varover,[],3);
segment_ind(:,:,2) = floor((vi-1)/ww); % vy
segment_ind(:,:,1) = mod((vi-1),ww);   % vx

[y, x] = meshgrid(1:1:width_cut, 1:1:height_cut);
x = x+ww-1;
y = y+ww-1;

segment_x_ind = zeros(height_cut,width_cut,wp);
segment_y_ind = zeros(height_cut,width_cut,wp);

for i=1:ww
    for j=1:ww
        segment_x_ind(:,:,i+(j-1)*ww) = x+segment_ind(:,:,1)-ww+i;
        segment_y_ind(:,:,i+(j-1)*ww) = y+segment_ind(:,:,2)-ww+j;
    end
end

ind_1d = (segment_y_ind(:)-1)*height + segment_x_ind(:);
im_1d = im(:);
segment_array_1d = im_1d(ind_1d);
segment_array = reshape(segment_array_1d,size(segment_y_ind));



if flat~=1
    gfilter_temp = zeros(1,1,wp);
    gfilter_temp(1,1,:) = gfilter(:);
    gfilter_ar = repmat(gfilter_temp,height_cut,width_cut,1);
    oim(ww:(height-ww+1),ww:(width-ww+1))=sum(segment.*gfilter_ar,3);
else
    oim(ww:(height-ww+1),ww:(width-ww+1))=mean(segment_array,3);
end


% for y=ww:(height-ww+1)
%     for x=ww:(width-ww+1)
%         [~, vi]=min(reshape(varnim((y-ws):(y+ws),(x-ws):(x+ws)),wp,1)); %the lowest variance v and vi its index
%         vi=vi-1;
%         vx=floor(vi/ww); %where it is inside the widow
%         vy=mod(vi,ww);
%         
%         segment=im((y+vy-ww+1):(y+vy),(x+vx-ww+1):(x+vx)); %cut the segment
%         if flat~=1
%             oim(y,x)=sum(sum(segment.*gfilter));
%         else
%             oim(y,x)=mean(mean(segment));
%         end
%     end
%     y
% end