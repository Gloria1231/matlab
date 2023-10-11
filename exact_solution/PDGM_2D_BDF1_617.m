function [phi,U,r0,ksi,psi0] = PDGM_2D_BDF1_617(N,T,dt)

global k2 kx ky L hx hy epsilon lamda D A  S1 S2 xi m tau K
            
%% 初始参数                    
M = round(T/dt);   % 时间划分
k=1;               % the indice of order

%区域
L=2*pi;            % domain
hx = L/N;          % space size
hy = hx;
x  =  hx*(0:N-1);  
y  =  x;

[xx,yy] = meshgrid(x,y); % mesh grid


phi0 = initphi(xx,yy);   % initial phase field
  U0  = initu(xx,yy);    % initial temperature field

[kx,ky,k2]=Fourier_2D(L,N);  % give the operator

energyflag = 0;          % the indice of plot energy
if 1 == energyflag
     filename_psi   = ['S1',num2str(S1),'S2',num2str(S2),'dt',num2str(dt),'_psi.txt'];
     filename_energy = ['S1',num2str(S1),'S2',num2str(S2),'dt',num2str(dt),'_energy.txt'];    
     out_psi = fopen(filename_psi,'w');
     out_energy = fopen(filename_energy,'w');
end

%初始时刻的r0
r0 =hx*hy*sum(sum(1/2*((fun_kappa(phi0).^2).*(-lap_diff(phi0).*phi0))+F(phi0) ...
    +lamda/(2*epsilon*K)*(U0.^2)+...
    +S2*1/(2*epsilon^2)*phi0.^2))+A;    % initial r0
ksi0=1;
psi0=1;
t=0;
%初始能量
if 1 == energyflag
    [energy1,energy2,err_energy,err_ksi] = calculate_energy(hx,hy,phi0,U0,r0,ksi0);
    fprintf(out_psi,'%14e  %18e %18e \n',t,psi0);
    fprintf(out_energy,'%14e  %18e  %18e  %18e  %18e \n',t,energy1,energy2,err_energy,err_ksi);
end

%% begin 
for nt = 1:M
    t = t+dt;
    
    phi_star = phi0;
    U_star   = U0;

    rhs1 = dt*fun_rhs1(xx,yy,t);                      
    rhs2 = dt*fun_rhs2(xx,yy,t);   

    % step 1
    g_phi    =tau*phi0-dt*S1*lap_diff(phi_star)+...         % tau*phin-dt*S1*\laplace phin
              +dt*var_der(phi_star)-dt*F_derivative(phi_star)+... % dt*\delta E_1/\delta \phi
              +dt*lamda/epsilon*P(phi_star).*U_star+rhs1;   % right hand term: lamda/epsilon*P*U=4*lamda*epsilon*F*U
    operator_phi=tau-S1*dt.*k2+dt*S2/(epsilon^2);           % left hand term：tau-S1*\laplace+S2/epsilon^2
    phi_hat  = real(ifft2(fft2(g_phi)./(operator_phi)));    % solve \phi_hat

    % step 2
    g_U     = U0-K*P(phi_star).*(phi_hat-phi0)+rhs2;        % right hand term: 
    U_hat   = real(ifft2(fft2(g_U)./(1-dt*D.*k2)));         % solve U_hat
     
    % Step 3   
    Energy_n = hx*hy*sum(sum(1/2*((fun_kappa(phi_hat).^2).*(-lap_diff(phi_hat).*phi_hat))+F(phi_hat) ...
               +lamda/(2*epsilon*K)*(U_hat.^2) +...
               +S2*1/(2*epsilon^2)*phi_hat.^2))+A;
    r_star = r0/(1+dt*fun_KK(phi_hat,U_hat)/Energy_n);       % compute q_hat
    
    % Step 4  
    ksi =  r_star/Energy_n;                                  % compute coefficient
    
    % Step 5
    eta = 1-(1-ksi)^(k+1);                                   % compute coefficient

    %% update phi, U
    phi = eta*phi_hat;
    U   = eta*U_hat; 
    

    % Step 6
    energy_new=hx*hy*sum(sum(1/2*((fun_kappa(phi).^2).*(-lap_diff(phi).*phi))+F(phi) ...
              +lamda/(2*epsilon*K)*(U.^2)+...
              +S2*1/(2*epsilon^2)*phi.^2))+A; 
    if r_star==energy_new 
        psi=0;
        gamma=r_star*ksi/(fun_KK(phi,U));
    elseif r_star>energy_new
        psi=0;
        gamma=(r_star-energy_new)/(dt*fun_KK(phi,U))+ksi*fun_KK(phi_hat,U_hat)/(fun_KK(phi,U));
    elseif r_star<energy_new && r_star-energy_new+dt*r_star*ksi*fun_KK(phi_hat,U_hat)>=0
        psi=0;
        gamma=(r_star-energy_new)/(dt*fun_KK(phi,U))+ksi*fun_KK(phi_hat,U_hat)/(fun_KK(phi,U));
    elseif r_star<energy_new && r_star-energy_new+dt*r_star*ksi*fun_KK(phi_hat,U_hat)<0
        psi=1-dt*ksi*fun_KK(phi_hat,U_hat)/(energy_new-r_star);
        gamma=0;
    end

    r=psi*r_star+(1-psi)*energy_new;                % Do relax
    r0=r;

    phi0=phi;
    U0=U;
      
    if 1 == energyflag
        [energy1,energy2,err_energy,err_ksi] = calculate_energy(hx,hy,phi,U,r,ksi);
         fprintf(out_psi,'%14e  %18e %18e \n',t,psi,eta);
         fprintf(out_energy,'%14e  %18e  %18e  %18e  %18e \n',t,energy1,energy2,err_energy,err_ksi);
    end  
end
if 1 == energyflag
    fclose(out_psi);
    fclose(out_energy);
end
end
