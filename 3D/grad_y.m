function grad_y=grad_y(phi)
global ky
grad_y=real(ifftn((ky.*fftn(phi))));
end