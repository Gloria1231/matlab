%%温度右端项
function z = fun_rhs2(x,y,t)
global K D
        z_exact=cos(x).*sin(y).*cos(t);
        z_phi_exact=sin(x).*cos(y).*cos(t);
        z_t=-cos(x).*sin(y).*sin(t);
        z_xx=-cos(x).*sin(y).*cos(t);
        z_phi_t=-sin(x).*cos(y).*sin(t);
        z_u=K*(1-z_phi_exact.^2).^2.*z_phi_t;
        z=z_t-2*D*z_xx-z_u;
 end