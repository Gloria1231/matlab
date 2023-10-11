function grad_y=grad_y(phi)
global ky
grad_y=real(ifft2((ky.*fft2(phi))));
end