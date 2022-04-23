function [K] = momentumMatrice_K_actin(parameter,cellinfo)

% VISCOUS MATRIX (Kv) -> -2*mu*int[dxWdxV]
% FRITION MATRIX (Kf) -> - eta*int[V*W]

%---------PARAMETERS FOR ISOPARAMETRIC INTEGRATION-------------------------
ngaus = cellinfo{1}.meshparam.ngaus ;  
p = cellinfo{1}.meshparam.p;  
[z, w] = quadrature_1D(ngaus);
[N,Nx,Nxx] = shapeFunc_1D(p,z);

%---------DEFINITION AND BUILDING OF STIFFNESS MATRICIES-------------------
X = cellinfo{1}.meshparam.X;
T = cellinfo{1}.meshparam.T;
nElem = cellinfo{1}.meshparam.Nele;
Kv = zeros(size(X,1),size(X,1)); %VISCOUS STIFFNESS MATRIX
Kf = zeros(size(X,1),size(X,1)); %FRICTION STIFFNESS MATRIX

% Loop on elements
for ielem = 1:nElem
    % Te: global number of the nodes in the current element
    Te = T(ielem,:);
    Xe = X(Te,:);
    adh_vec = [cellinfo{1}.meshparam.DOFadh_vec(Te(1,1),1)...
        cellinfo{1}.meshparam.DOFadh_vec(Te(1,2),1)];

    zeta_zero_vec = [parameter.zeta_m*(Te(1,1)) + parameter.zeta_b...
         parameter.zeta_m*(Te(1,2)) + parameter.zeta_b];
    % Element matrices, source terms and derivatives of residuals
    [Kfe,Kve] = EleMat_actin(Xe,cellinfo,parameter, N, Nx,w,adh_vec,zeta_zero_vec);
    % Assembly  

    r = [ielem ielem+1];
    Kv(r(1):r(1)+1,r(1):r(2)) = Kv(r(1):r(1)+1,r(1):r(2)) + Kve;
    Kf(r(1):r(1)+1,r(1):r(2)) = Kf(r(1):r(1)+1,r(1):r(2)) + Kfe;

end
K = Kv+Kf;

%------------DEFINITION OF ELEMENT MATRICIES-----------------------------
function  [Kfe,Kve] = EleMat_actin(Xe,cellinfo,parameter, N, Nx,w,adh_vec,zeta_zero_vec)
nen=cellinfo{1}.meshparam.nen;
ngaus = cellinfo{1}.meshparam.ngaus ; 
Kfe  = zeros(nen,nen);
Kve  = zeros(nen,nen);
visc = parameter.visc;
   
for ig = 1:ngaus
    N_ig    = N(ig,:);
    Nxi_ig  = Nx(ig,:);
    Jacob = Nxi_ig*Xe;
    dvolu = w(ig)*det(Jacob);
    res = Jacob\Nxi_ig;
    nx = res(1,:);
    adh_ig = N_ig*adh_vec';
    zeta_zero_ig = N_ig * zeta_zero_vec';
    
    zeta = parameter.zeta_small + (2*zeta_zero_ig*(adh_ig/parameter.Azeta)^2)/(1+(adh_ig/parameter.Azeta)^2);
    %zeta = parameter.zeta_small + (2*zeta_zero_ig)/(1+(adh_ig/parameter.Azeta)^2);
    
    Kfe = Kfe - zeta*(N_ig'*N_ig)*dvolu;
    Kve = Kve - visc*(nx'*nx)*dvolu;
end
