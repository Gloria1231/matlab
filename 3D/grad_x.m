function grad_x=grad_x(phi)
global kx
grad_x=real(ifftn((kx.*fftn(phi))));
end