%Solution of the equations for an AC grain (local problem)
%Explicit method (forward Euler)
% 1) Run the code and analyze the plots
% 2) Solve also the unloading problem: load for 200s, and unload (c_ext=0)
% for 100s

%BAX 1500
epse = 0.37; epsp = 0.8; rhoS = 2.22e6; K = 3.e-3; qm = 0.4; %Porosity constants, Laumuir's isoterm constants.
Dp=1.e-8; %Difusion intraparticular
B= 1.35;  %The constant related to the boundary condition (B.C) mass transfer. (small B, absorb slowly, big B, faster absortion)
R=0.002; %Radius of the bullet.

c_ext = 1300; %Exterior concentration

L=@(q) q./(K*(qm-q)); %Laumir's isoterm.
dL=@(q) qm./(K*(qm-q).^2);

Tfinal = 200;
dt = 0.01;  
nOfTimeSteps = round(Tfinal/dt);

alpha = Dp*B/R; beta=35*Dp/(R^2);
qb=zeros(nOfTimeSteps+1,1); qR=qb;

%Loop for time steps
for n=1:nOfTimeSteps
    cmL=c_ext-L(qR(n));
    qb(n+1)=qb(n)+dt*3*alpha*cmL;
    qR(n+1)=qR(n)+dt*(beta*(qb(n)-qR(n))+10*alpha*cmL);
end

ts = [0:nOfTimeSteps]*dt;
figure(1)
plot(ts,[qR,qb],'-*'), legend('qr','qb'),xlabel('t'),ylabel('q')

Tfinal = 100;
dt = 0.01;
nOfTimeStepsUnloading = round(Tfinal/dt);
qbLoading=qb; qRLoading=qR;
qb = zeros(nOfTimeStepsUnloading+1,1); qR = qb;
qb(1) = qbLoading(end); qR(1) = qRLoading(end);
c_ext = 0;

%Loop for time steps
for n=1:nOfTimeStepsUnloading
    cmL=c_ext-L(qR(n));
    qb(n+1)=qb(n)+dt*3*alpha*cmL;
    qR(n+1)=qR(n)+dt*(beta*(qb(n)-qR(n))+10*alpha*cmL);
end


ts = [0:nOfTimeStepsUnloading]*dt;
figure(2)
plot(ts,[qR,qb],'-*'), legend('qr','qb'),xlabel('t'),ylabel('q')

ts = [0:nOfTimeSteps + nOfTimeStepsUnloading+1]*dt;
figure(3)
qRT = [qRLoading;qR];
qbT = [qbLoading;qb];
plot(ts,[qRT,qbT], '-*'), legend('qr', 'qb'),xlabel('t'),ylabel('q')

% sigma = ( ((1-epse)/epse)*3*Dp*B/R )*((1-epsp)*rhoS-epsp*dL(qb));
% figure(2)
% plot(ts(2:end),sigma(2:end)), ylabel('sigma'), xlabel('t')
% 
