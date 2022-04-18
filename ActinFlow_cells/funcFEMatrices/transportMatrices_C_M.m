function  C = transportMatrices_C_M(cellinfo)

X = cellinfo{1}.meshparam.X;
T = cellinfo{1}.meshparam.T;
nElem = size(T,1); 
nen = cellinfo{1}.meshparam.nen;
nedof = nen; 
n = nedof^2*nElem; 
coef_C  = zeros(1,n);  indC_i  = zeros(1,n); indC_j  = zeros(1,n); indC = 1; 
% Loop on elements
for ielem = 1:nElem
    Te = T(ielem,:);
    Xe = X(Te,:);
    w_vel = [cellinfo{1}.meshparam.DOFw_vec(Te(1,1),1)...
        cellinfo{1}.meshparam.DOFw_vec(Te(1,2),1)];
    Ce = EleMat(Xe,cellinfo,w_vel);
    % Assembly
    Te_dof = Te; 
    for irow = 1:nedof
        for icol = 1:nedof
            indC_i(indC) = Te_dof(irow);
            indC_j(indC) = Te_dof(icol);
            coef_C(indC) = Ce(irow,icol);
            indC = indC+1;
        end
    end
end
indC_i = indC_i(1:indC-1);
indC_j = indC_j(1:indC-1);
coef_C = coef_C(1:indC-1);
C = sparse(indC_i,indC_j,coef_C);

function  Ce = EleMat(Xe,cellinfo,w_vel)
nen = cellinfo{1}.meshparam.nen;
Ce  = zeros(nen,nen);
p = cellinfo{1}.meshparam.p;                                         
ngaus = cellinfo{1}.meshparam.ngaus ;   
[z, w] = quadrature_1D(ngaus);
[N,Nx, Nxx] = shapeFunc_1D(p,z);

for ig = 1:ngaus
    N_ig    = N(ig,:);
    Nxi_ig  = Nx(ig,:);
    Jacob = Nxi_ig*Xe;
    dvolu = w(ig)*det(Jacob);
    res = Jacob\Nxi_ig;
    nx = res(1,:);
    
    v = N_ig*w_vel';
    Ce = Ce - v*nx'*N_ig*dvolu;
end



