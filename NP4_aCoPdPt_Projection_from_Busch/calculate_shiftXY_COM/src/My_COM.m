function [comX, comY] = My_COM(img)
    xdim = size(img,1);
    ydim = size(img,2);
    [Y, X] = meshgrid(1:1:ydim,1:1:xdim);
    comX = sum(sum(X.*img,2),1) / sum(img(:));
    comY = sum(sum(Y.*img,2),1) / sum(img(:));
end
    