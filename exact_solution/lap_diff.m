function lap=lap_diff(phi)
global k2
lap=real(ifft2((k2.*fft2(phi))));
end