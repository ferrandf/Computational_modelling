function errorL2=computeL2error(u,X,T,theReferenceElement)

IPweights = theReferenceElement.IPweights;
IPcoord = theReferenceElement.IPcoord;
N=theReferenceElement.N;
Nxi=theReferenceElement.Nxi;
Neta=theReferenceElement.Neta;

nOfElements = size(T,1);
errorL2 = 0;
%Loop in elements
for i=1:nOfElements
    Te=T(i,:); Xe=X(Te,:); 
    xe = Xe(:,1); ye=Xe(:,2);
    Ue=u(Te); 
    for g=1:length(IPweights)
        N_g = N(g,:);
        Nxi_g = Nxi(g,:);
        Neta_g = Neta(g,:);
        J = [Nxi_g*xe	  Nxi_g*ye
            Neta_g*xe  Neta_g*ye];
        dvolu=IPweights(g)*det(J);
        Xg = N_g*Xe;
        errorL2 = errorL2 + (analytical(Xg)-N_g*Ue)^2*dvolu;
    end      
end
errorL2 = sqrt(errorL2);
  