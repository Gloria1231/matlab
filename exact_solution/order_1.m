clc; clear
global  epsilon lamda D A S1 S2 xi m tau K
    epsilon =1;
    tau = 4e3;
    D = 1;
    lamda =1;    
    A = 1;
    K = 0.01; 
       
    xi = 0.05; 
    m= 4;
    S1=4;
    S2=4;
 
    N=128;
    %区域
    L=2*pi;            % domain
    hx = L/N;          % space size
    hy = hx;
    x  =  hx*(0:N-1);  
    y  =  x;

    [xx,yy] = meshgrid(x,y); % mesh grid
    T = 1;
    dt_arr = 0.2./2.^[1 2 3 4 5]'; 
    dt_leg=length(dt_arr);
    for k = 1:dt_leg
    dt=dt_arr(k);
    [phi,U,r0,ksi,psi0] = PDGM_2D_BDF1_617(N,T,dt);
    phi_exact = exactphi(xx,yy,T);
    U_exact = exactu(xx,yy,T);
    error1(k,1) = sqrt(sum(sum((phi_exact - phi).^2))*2*pi/N*2*pi/N);   % L2
    error2(k,1) = sqrt(sum(sum((U_exact - U).^2))*2*pi/N*2*pi/N);   % L2
    end

    %% Plot
    hh1=plot(log10(dt_arr(1:dt_leg)),log10(error1),'*-','LineWidth',1);
    hold on;
    hh2=plot(log10(dt_arr(1:dt_leg)),log10(error2),'o-','LineWidth',1);
    hold on;
    t=-2.2:.01:-1;
    yt=t+0.5;
    plot(t,yt,'r-.')
    xlabel('$log_{10}(\delta t)$','Interpreter','latex');
    ylabel('$log_{10}$(error)','Interpreter','latex');
    legend('$\phi$','$U$','$slope=1$','Interpreter','latex');
    title('BDF1','Interpreter','latex');
    grid on;
    hold on;