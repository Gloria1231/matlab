function grad_x=grad_x(phi)
global kx
grad_x=real(ifft2((kx.*fft2(phi))));
end