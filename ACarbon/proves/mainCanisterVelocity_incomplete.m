clear all
close all

%PREPROCES
%Reference element
degree = 1; typeOfElement=1; %1=TRI, 0=QUA
nOfElementNodes = 3;
theReferenceElement = createReferenceElement(degree,typeOfElement);
%Mesh: regular mesh in a rectangular domain [0,1]x[0,1]
nRef=1;
[X,T] = CreateMesh(typeOfElement,nOfElementNodes,[0,0.12,0,0.15],12*nRef+1,15*nRef+1);
[Xaux,Taux] = CreateMesh(typeOfElement,nOfElementNodes,[0.09,0.11,-0.02,0],2*nRef+1,2*nRef+1);
[X,T]=glueMeshes(X,T,Xaux,Taux);
[Xaux,Taux] = CreateMesh(typeOfElement,nOfElementNodes,[0.01,0.03,0.15,0.17],2*nRef+1,2*nRef+1);
[X,T]=glueMeshes(X,T,Xaux,Taux);
figure(2), clf
PlotMesh_lap(T,X,typeOfElement,'k-');

%Definition of Dirichlet boundary conditions (1st node set to 0)
x = X(:,1); y = X(:,2); tol=1.e-10;
nodesCCD = find((abs(y+0.02))<tol | abs(y-0.17)<tol);
nodestop = find(abs(X(:,2)-0.17)<tol);
nodesbot = find(abs(X(:,2)+0.02)<tol);
uCCD= [zeros(size(nodesbot));ones(size(nodestop))]; hold on, plot(x(nodesCCD),y(nodesCCD),'bo','MarkerSize',16); hold off
%Definition of connectivity matrices for non-homogeneous Neumann boundaries
Tboundary=connectivityMatrixBoundary(T,typeOfElement); %whole boundary
Xmid=(X(Tboundary(:,1),:)+X(Tboundary(:,2),:))/2;%mid points of the boundary sides
sidesTop=find(abs(Xmid(:,2)-0.17)<tol); 
TNeumannTop=Tboundary(sidesTop,:);
sidesBottom=find(abs(Xmid(:,2)+0.02)<tol);
TNeumannBottom=Tboundary(sidesBottom,:);
hold on, PlotBoundary([TNeumannTop;TNeumannBottom],X); hold off, axis off



%Definition of data functions
funcSourceTerm = @(X) zeros(size(X,1),1);
funcNeumann = @(X) ones(size(X,1),1);

[K,f]=computeSystemLaplace(X,T,theReferenceElement,funcSourceTerm);
%
%
%
% To Do: computation of potential and velocity at nodes ...
methodDBC = 1; % 1=System reduction or 2=Lagrange multipliers 
if methodDBC == 1 %System reduction
    unknowns= setdiff([1:size(X,1)],nodesCCD);
    f = f(unknowns)-K(unknowns,nodesCCD)*uCCD;
    K=K(unknowns,unknowns);
    %System solution
    sol=K\f;
    %Nodal values: system solution and Dirichlet values
    u = zeros(size(X,1),1);
    u(unknowns) = sol; u(nodesCCD) = uCCD;
else %LagrangeMultipliers
    nOfDirichletDOF = length(uCCD); nOfDOF = size(K,1);
    A = spalloc(nOfDirichletDOF,nOfDOF,nOfDirichletDOF);
    A(:,nodesCCD) = eye(nOfDirichletDOF);
    b = uCCD;
    Ktot = [K A'; A spalloc(nOfDirichletDOF,nOfDirichletDOF,0)];
    ftot = [f;b];
    sol = Ktot\ftot;
    u = sol(1:nOfDOF); lambda = sol(nOfDOF+1:end);
end

[ux,uy] = computeGradientSmoothing(u,X,T,theReferenceElement);
ux = ux/12;
uy = uy/12;
%
%

%streamLines
phi=computeStreamFunction(ux,uy,X,T,Tboundary,theReferenceElement);
figure(6), contourPlot(phi,X,T), title('Stream lines')
hold on, PlotBoundary(Tboundary,X); hold off, axis off
%magnitude of the velocity
figure(7), PlotNodalField(sqrt(ux.^2+uy.^2),X,T); title('|velocity|'), view(2), axis equal

%Saving the velocity at nodes
save velocityAtNodesNref1 ux uy X T typeOfElement degree 
