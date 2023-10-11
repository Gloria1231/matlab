%---------------------------------anisotropy coefficient--------------------------------
%fourfold symmetric model type anisotropy \kappa(grad(phi))
function kappa=fun_kappa(phi)
global xi m
angle=atan2(grad_y(phi),grad_x(phi));
kappa=1+xi*cos(m*angle);
end
