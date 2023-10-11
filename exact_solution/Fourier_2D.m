function [kx,ky,k2]= Fourier_2D(L,N)

%1i：复数单位
k1x = 1i*[0:N/2-1 0 -N/2+1:-1]*(2*pi/L);
k1y = k1x;

%以k1x为行，k1y为列的矩阵
[kx,  ky ] = meshgrid(k1x,k1y);

k2x = k1x.^2;
k2y = k1y.^2;
%定义laplace
[kxx, kyy] = meshgrid(k2x,k2y);

k2 = kxx + kyy;
end