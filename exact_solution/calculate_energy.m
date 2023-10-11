function [energy1,energy2,err_energy,err_ksi] = calculate_energy(hx,hy,phi,U,r,ksi)
global  epsilon lamda K A S1 S2
%原来的能量
energy1 = r-A;
energy2 = hx*hy*sum(sum( F(phi)))+S2*1/(2*epsilon^2)*hx*hy*sum(sum( phi.^2))+...
          0.5*hx*hy*sum(sum((fun_kappa(phi).^2).*(-lap_diff(phi).*phi)))+...
          lamda/(2*epsilon*K)* hx*hy*sum(sum(U.^2));

err_energy=abs(energy1-energy2);
err_ksi=abs(1-ksi);
end