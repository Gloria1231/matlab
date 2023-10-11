function kappa_a=fun_kappa_a(phi)
global xi m
angle=atan2(grad_y(phi),grad_x(phi));
kappa_a=-xi*sin(m*angle)*m;
end