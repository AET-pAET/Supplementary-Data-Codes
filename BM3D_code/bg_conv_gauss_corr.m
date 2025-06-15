% bg_conv_gauss_corr
% Author: Yongsoo Yang, UCLA Physics and Astronomy
%         yongsoo.ysyang@gmail.com

function after_conv = bg_conv_gauss_corr(bg_mask,convd,thresh)
    xdim = size(bg_mask,1);
    ydim = size(bg_mask,2);
    
    [Y, X] = meshgrid(0:1:ydim-1, 0:1:xdim-1);
    
    Kernel = exp(-1*((Y-ceil((ydim-1)/2)).^2+(X-ceil((xdim-1)/2)).^2)*log(2)/(convd)^2);
    
    mask_pad = zeros(2*xdim-1, 2*ydim-1);
    mask_pad(floor((xdim+1)/2):floor((xdim+1)/2)+xdim-1,floor((ydim+1)/2):floor((ydim+1)/2)+ydim-1) = bg_mask;
   
    Ker_pad = zeros(2*xdim-1, 2*ydim-1);
    Ker_pad(floor((xdim+1)/2):floor((xdim+1)/2)+xdim-1,floor((ydim+1)/2):floor((ydim+1)/2)+ydim-1) = Kernel;
    
    
    Conv = fftshift(ifftn(fftn(ifftshift(mask_pad)).*fftn(ifftshift(Ker_pad))));  
    
    after_conv = Conv(floor((xdim+1)/2):floor((xdim+1)/2)+xdim-1,floor((ydim+1)/2):floor((ydim+1)/2)+ydim-1);
    after_conv = after_conv > thresh;
end
                
                
