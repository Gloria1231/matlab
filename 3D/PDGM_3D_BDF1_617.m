function [phi,U,r0] = PDGM_3D_BDF1_617(N)
% Solve 3D phase field crystal equaiton
%    (1)    \tau \phi_t = \Delta \phi-f(phi)-lamda/epsilon*(1-\phi^2)^2*u
%    (2)    u_t = D*\Delta u + K*(1-\phi^2)^2*\phi_t
global k2 kx ky kz L hx hy hz epsilon lamda D A S1 S2 xi m tau K
  
%%
T  = 1800;          % total time 总时间
dt = 0.1;           % time step size 时间步长
M = round(T/dt);    % 时间划分

%参数
epsilon  = 3e-2;    % interface width
lamda = 260;        % linear kinetic coefficient
D     = 2e-4;       % diffusion rate
tau=2.5e4;          % mobility parameter
K=1.5;              % latent heat parameter
A=1;                % auxiliary variable constant

% 稳定化参数
S1=4;              % S1*\Laplace \phi: to enhace high order term
S2=4;              % S2*\phi         : to balance nonlinear term
xi=0.05;           % anisotropy strength
m=4;               % the number of anisotropy

k=1;               % the indice of order

% Space: Domain and N
N=128;             % Fourier modes
L = 2*pi;          % domain
hx = L/N;          % space size
hy = hx;
hz = hx;
x  =  hx*(0:N-1);
y  =  x;
z  =  x;

[xx,yy,zz] = meshgrid(x,y,z);   % mesh grid

%初值
phi0 = tanh((0.2-sqrt((xx-pi).^2+(yy-pi).^2+(zz-pi).^2))/0.072);  % initial phase field
  U0 = (phi0>0).*0+(phi0<=0).*(-0.55);                            % initial temperature field

% indice
    saveflag = 1;      
    energyflag = 1;  
% SAVE value
t = 0;
if 1 == saveflag
    if~isfolder('value')
    mkdir('value');
    end
    phi_sol = ['value' '/phi=' num2str(t) '.txt'];
    fid = fopen(phi_sol, 'wt');
    fprintf(fid, '%f\n', phi0(:));
    fclose(fid);
    U_sol = ['value' '/U=' num2str(t) '.txt'];
    fidd = fopen(U_sol, 'wt');
    fprintf(fidd, '%f\n', U0(:));
    fclose(fidd);
end

[kx,ky,kz,k2] = Fourier_3D(L,N);

if 1 == energyflag
    filename_psi   = ['e',num2str(epsilon),'_mass.txt'];
    filename_energy = ['e',num2str(epsilon),'_energy.txt'];    
    out_psi = fopen(filename_psi,'w');
    out_energy = fopen(filename_energy,'w');
end

%初始时刻的r0
r0 = hx*hy*hz*sum(sum(sum(1/2*((fun_kappa(phi0).^2).*(-lap_diff(phi0).*phi0))+F(phi0) ...
    +lamda/(2*epsilon*K)*(U0.^2)+...
    +S2*1/(2*epsilon^2)*phi0.^2)))+A;
ksi0=1;
psi0=0;
% Initial energy
if 1 == energyflag
    [energy1,energy2,err_energy,err_ksi] = calculate_energy(hx,hy,hz,phi0,U0,r0,ksi0);
    fprintf(out_psi,'%14e  %18e  \n',t,psi0);
    fprintf(out_energy,'%14e  %18e  %18e  %18e  %18e \n',t,energy1,energy2,err_energy,err_ksi);
end
for nt = 1:M
    t = t+dt;
    
    phi_star = phi0;
    U_star   = U0;
    
    % step 1
    g_phi    = tau*phi0-dt*S1*lap_diff(phi_star)+...
              +dt*var_der(phi_star)-dt*F_derivative(phi_star)+...
              -dt*lamda/epsilon*P(phi_star).*U_star;                        % tau*phin-dt*S1*\laplace phin
    phi_hat  = real(ifftn(fftn(g_phi)./(tau-S1*dt.*k2+dt*S2/(epsilon^2)))); % solve \phi_hat

    % step 2
    g_U     = U0+K*P(phi_star).*(phi_hat-phi0);                             % right hand term
    U_hat   = real(ifftn(fftn(g_U)./(1-dt*D.*k2)));                         % solve U_hat
    
    % Step 3 
    Energy_n = hx*hy*hz*sum(sum(sum(1/2*((fun_kappa(phi_hat).^2).*(-lap_diff(phi_hat).*phi_hat))+F(phi_hat) ...
               +lamda/(2*epsilon*K)*(U_hat.^2) +...
               +S2*1/(2*epsilon^2)*phi_hat.^2)))+A;
    r_star =  r0/(1+dt*fun_KK(phi_hat,U_hat)/Energy_n);                     % compute q_hat
    
    % Step 4  
    ksi =  r_star/Energy_n;                                                 % compute coefficient
    
    % Step 5
    eta = 1-(1-ksi)^(k+1);                                                  % compute coefficient

    %% update 
    phi = eta*phi_hat;
    U   = eta*U_hat; 

    % Step 6
    energy_new=hx*hy*hz*sum(sum(sum(1/2*((fun_kappa(phi).^2).*(-lap_diff(phi).*phi))+F(phi) ...
               +lamda/(2*epsilon*K)*(U.^2)+...
               +S2*1/(2*epsilon^2)*phi.^2)))+A;
    if r_star==energy_new 
        psi0=0;
        gamma=r_star*ksi/(fun_KK(phi,U));
    elseif r_star>energy_new
        psi0=0;
        gamma=(r_star-energy_new)/(dt*fun_KK(phi,U))+ksi*fun_KK(phi_hat,U_hat)/(fun_KK(phi,U));
    elseif r_star<energy_new && r_star-energy_new+dt*r_star*ksi*fun_KK(phi_hat,U_hat)>=0
        psi0=0;
        gamma=(r_star-energy_new)/(dt*fun_KK(phi,U))+ksi*fun_KK(phi_hat,U_hat)/(fun_KK(phi,U));
    elseif r_star<energy_new && r_star-energy_new+dt*r_star*ksi*fun_KK(phi_hat,U_hat)<0
        psi0=1-dt*ksi*fun_KK(phi_hat,U_hat)/(energy_new-r_star);
        gamma=0;
    end
    
    r=psi0*r_star+(1-psi0)*energy_new;                              % Do relax
    r0=r;

    phi0=phi;
    U0=U;
 
      
    if 1 == energyflag
        [energy1,energy2,err_energy,err_ksi] = calculate_energy(hx,hy,hz,phi,U,r,ksi);
         fprintf(out_psi,'%14e  %18e %18e \n',t,psi0,eta);
         fprintf(out_energy,'%14e  %18e  %18e  %18e  %18e \n',t,energy1,energy2,err_energy,err_ksi);
    end 
    if 0 == rem(nt,1/dt)
      if 1 == saveflag
            phi_sol = ['value' '/phi=' num2str(t) '.txt'];
            fid = fopen(phi_sol, 'wt');
            fprintf(fid, '%f\n', phi(:));
            fclose(fid);
            U_sol = ['value' '/U=' num2str(t) '.txt'];
            fidd = fopen(U_sol, 'wt');
            fprintf(fidd, '%f\n', U(:));
            fclose(fidd);
      end 
    end           
end

if 1 == energyflag
    fclose(out_psi);
    fclose(out_energy);
end