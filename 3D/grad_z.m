function grad_z=grad_z(phi)
global kz
grad_z=real(ifftn((kz.*fftn(phi))));
end