function fN = computefNeumannLinearApproximation(X,TN,normalFluxFunction)
%fN = computefNeumannLinearApproximation(X,TN,normalFluxFunction)
%Computes the r.h.s vector corresponding to the Neumann B.C.
% -n·grad(u)=normalFluxFunction

nOfNodes = size(X,1);
nOfElements = size(TN,1);

IPcoord = [-1/sqrt(3); 1/sqrt(3)]; 
IPweights = [1,1]; 
N = [(1 - IPcoord)/2   (1+IPcoord)/2]; 

fN=zeros(nOfNodes,1);
for i=1:nOfElements
    Te=TN(i,:); %index of the nodes in the element
    Xe=X(Te,:); %coordinates of the nodes in the element
    fNe=elementalComputations(normalFluxFunction,Xe,IPcoord,IPweights,N);
    fN(Te) = fN(Te) + fNe;
end

function fNe=elementalComputations(normalFluxFunction,Xe,IPcoord,IPweights,N)

nnodes = size(Xe,1);
fNe=zeros(nnodes,1);

x1 = Xe(1,:); x2 = Xe(2,:);
t = x2 - x1; h = norm(t);  
n = [t(2), -t(1)]; n = n/h; 
%Loop in integration points
for k=1:length(IPweights)
    Nk=N(k,:); %basis functions at integration point
    xk = Nk*Xe; %xy-coordinates of the integration point
    dl=IPweights(k)*h/2; %diferential of line
    fNe = fNe - Nk'*normalFluxFunction(xk)*dl;
end
  

