%FV solution of pure convection equation u_t + div(velocity u) = 0
clear all, clc, close all

%Problem data
finalTime=1000;
numericalFnfunction=@numericalNormalFluxPureConvection;
uInflow=0.1;

nOfComponents=1; smoothedPlots=1;

load velocityAtNodesNref1.mat;
%Faces information, etc
[elementsForSide,sidesForElement,Tsides,nOfInteriorSides,Ve,Ls,normalVectors,XmidVolumes,XmidSides] = FVpreprocess(X,T);
nOfVolumes=length(Ve); nOfSides=length(Ls); nOfExternalSides=nOfSides-nOfInteriorSides;
tol = 1.e-10;
nodestop = find(abs(X(:,2)-0.17)<tol);
nodesbot = find(abs(X(:,2)+0.02)<tol);

velocityAtNodes= (5.e-4/0.872032067541866)*[ux uy];%(5.e-4/0.872032067541866) *[ux uy];
velocitySides=(velocityAtNodes(Tsides(:,1),:)+velocityAtNodes(Tsides(:,2),:))/2;
velocityDotNormals=velocitySides(:,1).*normalVectors(1,:)'+velocitySides(:,2).*normalVectors(2,:)';

%Identification of inflow and outflow boundary
tol=1.e-10; flowSides=find(abs(XmidSides(:,2)+0.02)<tol | abs(XmidSides(:,2)-0.17)<tol )';
impermeableSides=setdiff(nOfInteriorSides+[1:nOfExternalSides],flowSides);
figure(1), hold on, plot(XmidSides(flowSides,1),XmidSides(flowSides,2),'r*',XmidSides(impermeableSides,1),XmidSides(impermeableSides,2),'b*'), hold off

%Initial Condition
u=zeros(nOfComponents,nOfVolumes);

mindx=sqrt(min(Ve)*2);
dt=0.4*mindx/max(sqrt(velocitySides(:,1).^2+velocitySides(:,2).^2)); %time step, for Courant(=velocity*dt/dx)<=1


nOfTimeSteps=round(finalTime/dt);
U=u(1:nOfVolumes); nStepsU=1;%round(nOfTimeSteps/10);
nCompVe=repmat(Ve,nOfComponents,1);
%__________________________________________________________________________
%Loop in time steps
for n=1:nOfTimeSteps
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
    u = u - dt*numericalNormalFluxes./nCompVe;
    if mod(n,nStepsU)==0, U=[U;u(1,:)]; end
    u(u>uInflow)=uInflow; u(u<0)=0;
end

%__________________________________________________________________________
%Plots of the solution
if smoothedPlots==0
    for k=2:size(U,1)
        plotFVsolution(X,T,U(k,:),2), axis equal %, view(2), colorbar
        title(sprintf('t=%g',(k-1)*nStepsU*dt))
        kStep=(k-1)*nStepsU;
        if kStep<10, print(sprintf('plots/convection0%d',kStep),'-djpeg'), else, print(sprintf('plots/convection%d',kStep),'-djpeg'), end
    end
else
    for k=2:size(U,1)
        uNodes=computeNodesMeanValue(U(k,:),T);
        trisurf(T,X(:,1),X(:,2),uNodes), axis equal, shading interp %, view(2), colorbar
        title(sprintf('t=%g',(k-1)*nStepsU*dt))
        kStep=(k-1)*nStepsU;
        if kStep<10, print(sprintf('plots/convection0%d',kStep),'-djpeg'), else, print(sprintf('plots/convection%d',kStep),'-djpeg'), end
    end
end
view(2)
