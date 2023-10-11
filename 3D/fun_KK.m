function KK = fun_KK(phi,U)
global hx hy hz tau lamda D epsilon K S1 S2
E1 = 1/tau*((var_der(phi)-F_derivative(phi)+ ...
        - S2/(epsilon^2)*phi-lamda/epsilon*P(phi).*U).^2);
E2 = lamda*D/(epsilon*K)*(-lap_diff(U).*U);
KK  = hx*hy*hz*sum(sum(sum(E1+E2)));
end