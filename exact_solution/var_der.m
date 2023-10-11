%H(phi)=delta\kappa(grad_phi)/delta(phi)
function var_der=var_der(phi)
global S1
var_der=grad_y(fun_kappa(phi).*fun_kappa_a(phi).*grad_x(phi))+...
-grad_x(fun_kappa(phi).*fun_kappa_a(phi).*grad_y(phi))+...
+grad_x(fun_kappa(phi).^2.*grad_x(phi))+grad_y(fun_kappa(phi).^2.*grad_y(phi));
end