function  [M,K] = transportMatrices_MK_signal(cellinfo)

X = cellinfo{1}.meshparam.X;
T = cellinfo{1}.meshparam.T;
nElem = cellinfo{1}.meshparam.Nele;
nen=cellinfo{1}.meshparam.nen;
nedof = nen;  
n = nedof^2*nElem; 
coef_K  = zeros(1,n);  indK_i  = zeros(1,n); indK_j  = zeros(1,n); indK = 1; 
coef_M  = zeros(1,n);  indM_i  = zeros(1,n); indM_j  = zeros(1,n); indM = 1; 

% Loop on elements
for ielem = 1:nElem
    % Te: global number of the nodes in the current element
    Te = T(ielem,:);
    Xe = X(Te,:);
        
    % Element matrices, source terms and derivatives of residuals
    [Me,Ke] = EleMat(Xe,cellinfo);
    % Assembly
    Te_dof = Te; 
    for irow = 1:nedof
        for icol = 1:nedof
            indK_i(indK) = Te_dof(irow);
            indK_j(indK) = Te_dof(icol);
            coef_K(indK) = Ke(irow,icol);
            indK = indK+1;
            
            indM_i(indM) = Te_dof(irow);
            indM_j(indM) = Te_dof(icol);
            coef_M(indM) = Me(irow,icol);
            indM = indM+1;

        end
    end
end
indK_i = indK_i(1:indK-1);
indK_j = indK_j(1:indK-1);
coef_K = coef_K(1:indK-1);
K = sparse(indK_i,indK_j,coef_K);

indM_i = indM_i(1:indM-1);
indM_j = indM_j(1:indM-1);
coef_M = coef_M(1:indM-1);
M = sparse(indM_i,indM_j,coef_M);






function  [Me,Ke] = EleMat(Xe,cellinfo)
nen=cellinfo{1}.meshparam.nen;
Me  = zeros(nen,nen);
Ke  = zeros(nen,nen);

ngaus = cellinfo{1}.meshparam.ngaus ;  
p = cellinfo{1}.meshparam.p;  
[z, w] = quadrature_1D(ngaus);
[N,Nx,Nxx] = shapeFunc_1D(p,z);

for ig = 1:ngaus
    N_ig    = N(ig,:);
    Nxi_ig  = Nx(ig,:);
    Jacob = Nxi_ig*Xe;
    dvolu = w(ig)*det(Jacob);
    res = Jacob\Nxi_ig;
    nx = res(1,:);

    Me = Me + N_ig'*N_ig*dvolu; 
    Ke = Ke - (nx'*nx)*dvolu; 
end



