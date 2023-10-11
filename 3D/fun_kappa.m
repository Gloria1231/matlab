function kappa=fun_kappa(phi)
global xi 
grad_4=(grad_x(phi).^2+grad_y(phi).^2+grad_z(phi).^2).^2+1e-9;
r1=(grad_x(phi).^4+grad_y(phi).^4+grad_z(phi).^4)./grad_4;
a_1=(1-3*xi);
b_1=4*xi/(1-3*xi);
kappa=a_1*(1+4*xi/(1-3*xi)*r1);
end