%FV solution of pure convection equation u_t + div(velocity u) = 0
%clear all, clc, close all

%Problem data
finalTime=100;
velocityFunction=@(X) 15*[ux uy]; %definition of the velocity
numericalFnfunction=@numericalNormalFluxPureConvection;
uInflow=1;

nOfComponents=1; smoothedPlots=1;

%Mesh (X: nodal coordinates, T: connectivity matrix defining the elements/volumes)
nRef=1; 
%[X,T] = createRectangleMesh(1,3,0,10,0,0.5,20*nRef,nRef);
typeOfElement = 1;
nOfElementNodes = 3;
nOfElem1d = 20;
%[X,T] = createRectangleMesh(typeOfElement,nOfElementNodes,0,12,0,15,nOfElem1d,nOfElem1d); %Definition of the mesh

pgon = polyshape([0 7.8 7.8 8.2 8.2 12 12 4.2 4.2 3.8 3.8 0],[0 0 10 10 0 0 15 15 5 5 15 15]);
tr = triangulation(pgon);
tnodes = tr.Points';
telements = tr.ConnectivityList';
model = createpde;
geometryFromMesh(model,tnodes,telements);
pdegplot(model)
m1 = generateMesh(model, 'Hmax', 0.4, 'GeometricOrder','linear');
X = m1.Nodes; X = X';
T = m1.Elements; T = T';



figure(1), plotMesh(X,T)
%Faces information, etc
[elementsForSide,sidesForElement,Tsides,nOfInteriorSides,Ve,Ls,normalVectors,XmidVolumes,XmidSides] = FVpreprocess(X,T);
nOfVolumes=length(Ve); nOfSides=length(Ls); nOfExternalSides=nOfSides-nOfInteriorSides;

velocityAtNodes=velocityFunction(X);
velocitySides=(velocityAtNodes(Tsides(:,1),:)+velocityAtNodes(Tsides(:,2),:))/2;
velocityDotNormals=velocitySides(:,1).*normalVectors(1,:)'+velocitySides(:,2).*normalVectors(2,:)';

%Identification of inflow and outflow boundary

tol=1.e-10; flowSides=find( (abs(XmidSides(:,2))<tol & abs(XmidSides(:,1)-10.25) < 1)| (abs(XmidSides(:,2)-15)<tol & abs(XmidSides(:,1)-2) < 1))';
impermeableSides=setdiff(nOfInteriorSides+[1:nOfExternalSides],flowSides);
figure(1), hold on, plot(XmidSides(flowSides,1),XmidSides(flowSides,2),'r*',XmidSides(impermeableSides,1),XmidSides(impermeableSides,2),'b*'), hold off

%Initial Condition
u=zeros(nOfComponents,nOfVolumes);

mindx=sqrt(min(Ve)*2);
dt=1.3*mindx/max(sqrt(velocitySides(:,1).^2+velocitySides(:,2).^2)); %time step, for Courant(=velocity*dt/dx)<=1

%Try also dt=0.5*...  and dt=0.6*...

nOfTimeSteps=round(finalTime/dt);
U=u(1:nOfVolumes); nStepsU=1;%round(nOfTimeSteps/10);
nCompVe=repmat(Ve,nOfComponents,1);
%__________________________________________________________________________
%Loop in time steps
for n=1:nOfTimeSteps
    numericalNormalFluxes=zeros(nOfComponents,nOfVolumes,1);
    %Loop in internal sides/faces
    for s=1:nOfInteriorSides
        eL=elementsForSide(s,1); eR=elementsForSide(s,2); %The two elements that share the side
        Fn=numericalFnfunction(u(eL),u(eR),velocityDotNormals(s)); %Compute the flux using the values on elements L and R.
        numericalNormalFluxes(:,eL)=numericalNormalFluxes(:,eL)+Ls(s)*Fn; %Ls is the length of the side
        numericalNormalFluxes(:,eR)=numericalNormalFluxes(:,eR)-Ls(s)*Fn;
    end
    %Loop in inflow and outflow sides (flux=0 on impermeable sides) The
    %boundary conditions fix this values
    for s=flowSides
        eL=elementsForSide(s,1);
        Fn=numericalFnfunction(u(eL),uInflow,velocityDotNormals(s)); %The uInflow is the value set on the boundary.
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
    for k=2:20:size(U,1)
        plotFVsolution(X,T,U(k,:),2), axis equal %, view(2), colorbar
        title(sprintf('t=%g',(k-1)*nStepsU*dt))
        kStep=(k-1)*nStepsU;
        if kStep<10, print(sprintf('plots/convection0%d',kStep),'-djpeg'), else, print(sprintf('plots/convection%d',kStep),'-djpeg'), end
    end
else
    for k=2:20:size(U,1)
        uNodes=computeNodesMeanValue(U(k,:),T);
        trisurf(T,X(:,1),X(:,2),uNodes), axis equal, shading interp %, view(2), colorbar
        title(sprintf('t=%g',(k-1)*nStepsU*dt))
        kStep=(k-1)*nStepsU;
        if kStep<10, print(sprintf('plots/convection0%d',kStep),'-djpeg'), else, print(sprintf('plots/convection%d',kStep),'-djpeg'), end
    end
end
view(2)
