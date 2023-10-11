function [kx,ky,kz,k2] = Fourier_3D(L,N)

k1x = 1i*[0:N/2-1 0 -N/2+1:-1]*(2*pi/L);
k1y = k1x;
k1z = k1x;

[kx, ky, kz] = meshgrid(k1x,k1y,k1z);

k2x = k1x.^2;
k2y = k1y.^2;
k2z = k1z.^2;

[kxx, kyy, kzz] = meshgrid(k2x,k2y,k2z);

k2 = kxx + kyy + kzz;
end 