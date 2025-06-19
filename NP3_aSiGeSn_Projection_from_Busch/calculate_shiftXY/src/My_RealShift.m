function img2 = My_RealShift(img,dx,dy)

ny = size(img,1); nx = size(img,2);
[X, Y] = meshgrid(1:nx,1:ny);
X1=mod(X-dy-1,nx)+1;
Y1=mod(Y-dx-1,ny)+1;

img(:,end+1)=img(:,1);
img(end+1,:)=img(1,:);
ny = size(img,1); nx = size(img,2);
[X, Y] = meshgrid(1:nx,1:ny);
img2 = interp2(X,Y,img,X1,Y1,'cubic',0);

end