%Simulation of an AC filter, 1D problem. FD solution with upwind for
%convective term and forward Euler in time. It solves both global and local
%problem.
close all

%Material parameters
epse = 0.37; epsp = 0.8; 
rhoS = 2.22e6; 
K = 3.e-3; qm = 0.4;
Dp=1.e-8;
Kf = 0.6e-2; %0.6 cm/s
B=Kf/(Dp*rhoS*(1-epsp));
R=0.002;

c_ext = 1300;
velocity = 5.e-4;%cm/s
nu=1.e-9;
canisterSize=0.1;%10 cm

L=@(q) q./(K*(qm-q));
dL=@(q) qm./(K*(qm-q).^2);
sigmaConstant = ((1-epse)/epse)*3*Dp*B/R;
alpha = Dp*B/R; beta=35*Dp/(R^2);

N=50; %number of segments
dx = canisterSize/N;
c=zeros(N+2,1); c(1)=c_ext; %including a fictitius node for the right boundary (homogeneous Neumann)

%Tfinal=100; dt=2; fprintf('\nCourant number=%0.2f\n',velocity*dt/dx); %Final time and step for convection-diffusion
Tfinal=30000; dt=0.01; %Final time and time step for AC filter

nOfTimeSteps=round(Tfinal/dt);
C=zeros(N+1,nOfTimeSteps+1); Qb=C; QR=C;
C(:,1)=c(1:end-1);

for n=1:nOfTimeSteps
    i=2:N+1;
    sigma = sigmaConstant*((1-epsp)*rhoS+epsp*dL(Qb(i,n)));
    c(i)=c(i)+dt*(velocity*(c(i-1)-c(i))/dx+nu*(c(i+1)-2*c(i)+c(i-1))/(dx^2))+dt*sigma.*(L(QR(i,n))-c(i));
    c(1)=c_ext; c(end)=c(end-2);
    C(:,n+1)=c(1:end-1);
    cmL=C(:,n)-L(QR(:,n));
    Qb(:,n+1)=Qb(:,n)+dt*3*alpha*cmL;
    QR(:,n+1)=QR(:,n)+dt*(beta*(Qb(:,n)-QR(:,n))+10*alpha*cmL);
end
x=0:dx:canisterSize;
tplots=1:round(nOfTimeSteps/10):nOfTimeSteps+1;
figure(1), plot(x,C(:,tplots)); title("dt = " + dt), xlabel('x'), ylabel('c')
figure(2), plot(x,Qb(:,tplots)); title('qb')
figure(3), plot(x,QR(:,tplots)); title('qR')

Cl=C; Qbl=Qb; QRl=QR;

c_ext = 0;
velocity = -5.e-2;%cm/s
nu=1.e-9;
canisterSize=0.1;%10 cm

L=@(q) q./(K*(qm-q));
dL=@(q) qm./(K*(qm-q).^2);
sigmaConstant = ((1-epse)/epse)*3*Dp*B/R;
alpha = Dp*B/R; beta=35*Dp/(R^2);

N=50; %number of segments
dx = canisterSize/N;
c = zeros(N+2,1);
c = [Cl(2,end);Cl(:,end)]; c(N+2)=c_ext; %including a fictitius node for the left boundary (homogeneous Neumann)

%Tfinal=100; dt=2; fprintf('\nCourant number=%0.2f\n',velocity*dt/dx); %Final time and step for convection-diffusion
Tfinal=2000; dt=0.001;%0.033; %Final time and time step for AC filter

nOfTimeSteps=round(Tfinal/dt);
C=zeros(N+1,nOfTimeSteps+1); Qb=C; QR=C;
C(:,1)=c(2:end); Qb(:,1) = Qbl(:,end); QR(:,1) = QRl(:,end);

for n=1:nOfTimeSteps
    i=N+1:-1:2;
    sigma = sigmaConstant*((1-epsp)*rhoS+epsp*dL(Qb(i,n)));
    c(i)=c(i)+dt*(velocity*(c(i)-c(i+1))/dx+nu*(c(i-1)-2*c(i)+c(i+1))/(dx^2))+dt*sigma.*(L(QR(i,n))-c(i));
    c(end)=c_ext; c(1)=c(3);
    C(:,n+1)=c(2:end);
    cmL=C(:,n)-L(QR(:,n));
    Qb(:,n+1)=Qb(:,n)+dt*3*alpha*cmL;
    QR(:,n+1)=QR(:,n)+dt*(beta*(Qb(:,n)-QR(:,n))+10*alpha*cmL);
end
x=0:dx:canisterSize;
tplots=1:round(nOfTimeSteps/10):nOfTimeSteps+1;
figure(4), plot(x,C(:,tplots)); title("dt = " + dt), xlabel('x'), ylabel('c')
figure(5), plot(x,Qb(:,tplots)); title('qb')
figure(6), plot(x,QR(:,tplots)); title('qR')
