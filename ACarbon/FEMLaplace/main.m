%PREPROCES
%Reference element
degree = 1; typeOfElement=0; %1=TRI, 0=QUA
theReferenceElement = createReferenceElement(degree,typeOfElement);
nOfElementNodes = size(theReferenceElement.N,2);
figure(1), drawReferenceElement(theReferenceElement);
%Mesh: regular mesh in a rectangular domain [0,1]x[0,1]
nOfElem1d=4; %20
[X,T] = CreateMesh(typeOfElement,nOfElementNodes,[0,1,0,1],nOfElem1d,nOfElem1d);
figure(2), clf
PlotMesh(T,X,typeOfElement,'k-',1);
%Definition of Dirichlet boundary conditions on {x=0}U{y=0}U{y=1}
x = X(:,1); y = X(:,2); tol=1.e-10;
nodesCCD = find(abs(x)<tol|abs(y)<tol|abs(y-1)<tol); %Nodes on the boundary
hold on, plot(x(nodesCCD),y(nodesCCD),'bo','MarkerSize',16); hold off
uCCD=DirichletValue(X(nodesCCD,:));
%Definition of connectivity matrix for Neumann boundary {x=1}
nodesNeumann = find(abs(x-1)<tol); 
[kk,permut]=sort(X(nodesNeumann,2)); %ordering with increasing y
nodesNeumann=nodesNeumann(permut);
TNeumann=[nodesNeumann(1:end-1),nodesNeumann(2:end)];

%System of equations (without BC)
[K,f]=computeSystemLaplace(X,T,theReferenceElement,@sourceTerm);

%Neumann boundary conditions
f = f + computefNeumannLinearApproximation(X,TNeumann,@NeumannFunction);

%Dirichlet Boundary conditions & SYSTEM SOLUTION
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

%POSTPROCESS
figure(3)
PlotNodalField(u,X,T), title('FEM solution')
[Xpg,Fgp] = gradientElementalField(u,X,T,theReferenceElement);
figure(4)
PlotMesh(T,X,typeOfElement,'k-'); hold on, quiver(Xpg(:,1),Xpg(:,2),Fgp(:,1),Fgp(:,1),'LineWidth',2), hold off
%Gradient smoothing (L2 projection of the gradient onto de FE space)
[ux,uy] = computeGradientSmoothing(u,X,T,theReferenceElement);
figure(5)
PlotMesh(T,X,typeOfElement,'k-'); hold on, quiver(X(:,1),X(:,2),ux,uy,'LineWidth',2), hold off


%Comparison with analytical solution (if available)
figure(11)
[Xfine,Tfine] = CreateMesh(typeOfElement,nOfElementNodes,[0,1,0,1],40,40);
PlotNodalField(analytical(Xfine),Xfine,Tfine), title('Analytical solution')
L2error=computeL2error(u,X,T,theReferenceElement)


