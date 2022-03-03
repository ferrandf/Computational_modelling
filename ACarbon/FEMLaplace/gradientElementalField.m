function [Xpg,Fgp] = gradientElementalField(u,X,T,theReferenceElement)
% [Xpg,Fgp] = gradientElementalField(u,X,T,theReferenceElement)
% Computes the gradient of u at Gauss points from nodal values of u
% Xpg:  coodinates of Gauss points
% Fgp:  gradient at Gauss points
% quiver(Xpg(:,1),Xpg(:,2),Fgp(:,1),Fgp(:,1)) plots the field


IPweights = theReferenceElement.IPweights;
IPcoord = theReferenceElement.IPcoord;
N=theReferenceElement.N;
Nxi=theReferenceElement.Nxi;
Neta=theReferenceElement.Neta;

nOfElements  = size(T,1);
nOfIPElem = size(IPcoord,1);
nOfIP = nOfElements*nOfIPElem;

Xpg = zeros(nOfIP,2);
Fgp = zeros(nOfIP,2);

%Loop in elements
ipg=1;
for ielem = 1:nOfElements
  Te = T(ielem,:);  Xe = X(Te,:);
  %Loop in Gauss points
  for igaus = 1:nOfIPElem
    N_igaus = N(igaus,:); 
    Nxi_igaus = Nxi(igaus,:);  
    Neta_igaus = Neta(igaus,:);  
    Jacob = [Nxi_igaus*(Xe(:,1))	Nxi_igaus*(Xe(:,2))  
             Neta_igaus*(Xe(:,1))	Neta_igaus*(Xe(:,2))];
    T_xi  = Nxi_igaus*u(Te);
    T_eta = Neta_igaus*u(Te);
    res = Jacob\[T_xi;T_eta];
    u_x = res(1); u_y = res(2);
    Fgp(ipg,:) = [u_x,u_y];
    Xpg(ipg,:) = N_igaus*Xe;
    ipg = ipg + 1;
  end
end
