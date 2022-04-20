%FV solution of pure convection equation u_t + div(velocity u) = 0
clear all, clc, close all

%Problem data
finalTime=1000;
numericalFnfunction=@numericalNormalFluxPureConvection;
uInflow=0.1;

%Material parameters
epse = 0.37; epsp = 0.8; 
rhoS = 2.22e6; 
K = 3.e-3; qm = 0.4;
Dp=1.e-8;
Kf = 0.6e-2; %0.6 cm/s
B=Kf/(Dp*rhoS*(1-epsp));
R=0.002;
nu=1.e-9;


L=@(q) q./(K*(qm-q));
dL=@(q) qm./(K*(qm-q).^2);
sigmaConstant = 0;%((1-epse)/epse)*3*Dp*B/R;
alpha = Dp*B/R; beta=35*Dp/(R^2);

nOfComponents=1; smoothedPlots=0;

load velocityAtNodesNref1.mat;
%Faces information, etc
[elementsForSide,sidesForElement,Tsides,nOfInteriorSides,Ve,Ls,normalVectors,XmidVolumes,XmidSides] = FVpreprocess(X,T);
nOfVolumes=length(Ve); nOfSides=length(Ls); nOfExternalSides=nOfSides-nOfInteriorSides;
tol = 1.e-10;
nodestop = find(abs(X(:,2)-0.17)<tol);
nodesbot = find(abs(X(:,2)+0.02)<tol);

velocityAtNodes= (5.e-4/0.97)*[ux uy];%(5.e-4/0.872032067541866) *[ux uy];
velocitySides=(velocityAtNodes(Tsides(:,1),:)+velocityAtNodes(Tsides(:,2),:))/2;
velocityDotNormals=velocitySides(:,1).*normalVectors(1,:)'+velocitySides(:,2).*normalVectors(2,:)';

%Identification of inflow and outflow boundary
tol=1.e-10; flowSides=find(abs(XmidSides(:,2)+0.02)<tol | abs(XmidSides(:,2)-0.17)<tol )';
impermeableSides=setdiff(nOfInteriorSides+[1:nOfExternalSides],flowSides);
figure(1), hold on, plot(XmidSides(flowSides,1),XmidSides(flowSides,2),'r*',XmidSides(impermeableSides,1),XmidSides(impermeableSides,2),'b*'), hold off

%Initial Condition
u=zeros(nOfComponents,nOfVolumes);

mindx=sqrt(min(Ve)*2);
dt=0.0065*mindx/max(sqrt(velocitySides(:,1).^2+velocitySides(:,2).^2)) %time step, for Courant(=velocity*dt/dx)<=1


nOfTimeSteps=round(finalTime/dt)
U=u(1:nOfVolumes); nStepsU=1;%round(nOfTimeSteps/10);
nCompVe=repmat(Ve,nOfComponents,1);
qb=zeros(nOfVolumes,1); qr=qb;

%__________________________________________________________________________
%Loop in time steps
for n=1:nOfTimeSteps
    if mod(n,500)==0
        nOfTimeSteps - n
    end
    sigma = sigmaConstant*((1-epsp)*rhoS+epsp*dL(qb));
    numericalNormalFluxes=zeros(nOfComponents,nOfVolumes,1);
    %Loop in internal sides
    for s=1:nOfInteriorSides
        eL=elementsForSide(s,1); eR=elementsForSide(s,2);
        Fn=numericalFnfunction(u(eL),u(eR),velocityDotNormals(s));
        numericalNormalFluxes(:,eL)=numericalNormalFluxes(:,eL)+Ls(s)*Fn;
        numericalNormalFluxes(:,eR)=numericalNormalFluxes(:,eR)-Ls(s)*Fn;
    end
    %Loop in inflow and outflow sides (flux=0 on impermeable sides)
    for s=flowSides
        eL=elementsForSide(s,1);
        Fn=numericalFnfunction(u(eL),uInflow,velocityDotNormals(s));
        numericalNormalFluxes(:,eL)=numericalNormalFluxes(:,eL)+Ls(s)*Fn;
    end
    %Update of value at volumes
    u = u - dt*numericalNormalFluxes./nCompVe +dt*sigma.*(L(qr)-u);
    if mod(n,nStepsU*10)==0, U=[U;u(1,:)]; end
    u(u>uInflow)=uInflow; u(u<0)=0;
    cmL=u-L(qr);
    qr=qr+dt*(beta*(qb-qr)+10*alpha*cmL);
    qb=qb+dt*3*alpha*cmL;
    
    
    
end

% 
% 
% 
% 
% for n=1:nOfTimeSteps
%     i=2:N+1;
%     sigma = sigmaConstant*((1-epsp)*rhoS+epsp*dL(Qb(i,n)));
%     c(i)=c(i)+dt*(velocity*(c(i-1)-c(i))/dx+nu*(c(i-1)-2*c(i)+c(i+1))/(dx^2))+dt*sigma.*(L(QR(i,n))-c(i));
%     c(1)=c_ext; c(end)=c(end-2);
%     C(:,n+1)=c(1:end-1);
%     cmL=C(:,n)-L(QR(:,n));
%     Qb(:,n+1)=Qb(:,n)+dt*3*alpha*cmL;
%     QR(:,n+1)=QR(:,n)+dt*(beta*(Qb(:,n)-QR(:,n))+10*alpha*cmL);
% end
% x=0:dx:canisterSize;
% tplots=1:round(nOfTimeSteps/10):nOfTimeSteps+1;
% figure(1), plot(x,C(:,tplots)); title("Dp = " + Dp + " dx = " + dx + " dt = " + dt), xlabel('x'), ylabel('c')
% figure(2), plot(x,Qb(:,tplots)); title('qb')
% figure(3), plot(x,QR(:,tplots)); title('qR')
% 
% Cl=C; Qbl=Qb; QRl=QR;



%__________________________________________________________________________
%Plots of the solution

if smoothedPlots==0
    for k=2:60:size(U,1)
        plotFVsolution(X,T,U(k,:),2), axis equal %, view(2), colorbar
        title(sprintf('t=%g',(k-1)*nStepsU*dt))
        kStep=(k-1)*nStepsU;
        if kStep<10, print(sprintf('plots/convection0%d',kStep),'-djpeg'), else, print(sprintf('plots/convection%d',kStep),'-djpeg'), end
    end
else
    for k=2:60:size(U,1)
        uNodes=computeNodesMeanValue(U(k,:),T);
        trisurf(T,X(:,1),X(:,2),uNodes), axis equal, shading interp %, view(2), colorbar
        title(sprintf('t=%g',(k-1)*nStepsU*dt))
        kStep=(k-1)*nStepsU;
        if kStep<10, print(sprintf('plots/convection0%d',kStep),'-djpeg'), else, print(sprintf('plots/convection%d',kStep),'-djpeg'), end
    end
end
view(2)
