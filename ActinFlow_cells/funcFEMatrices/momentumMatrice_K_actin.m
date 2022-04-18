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
    % Element matrices, source terms and derivatives of residuals
    [Kfe,Kve] = EleMat_actin(Xe,cellinfo,parameter, N, Nx,w);
    % Assembly  

    r = [ielem ielem+1];
    Kv(r(1):r(1)+1,r(1):r(2)) = Kv(r(1):r(1)+1,r(1):r(2)) + Kve;
    Kf(r(1):r(1)+1,r(1):r(2)) = Kf(r(1):r(1)+1,r(1):r(2)) + Kfe;

end
K = Kv+Kf;

%------------DEFINITION OF ELEMENT MATRICIES-----------------------------
function  [Kfe,Kve] = EleMat_actin(Xe,cellinfo,parameter, N, Nx,w)
nen=cellinfo{1}.meshparam.nen;
ngaus = cellinfo{1}.meshparam.ngaus ; 
Kfe  = zeros(nen,nen);
Kve  = zeros(nen,nen);
visc = parameter.visc;
eta = parameter.eta;
eta_x = eta*ones(ngaus,1);
   
for ig = 1:ngaus
    N_ig    = N(ig,:);
    Nxi_ig  = Nx(ig,:);
    Jacob = Nxi_ig*Xe;
    dvolu = w(ig)*det(Jacob);
    res = Jacob\Nxi_ig;
    nx = res(1,:);
   
    
    Kfe = Kfe - eta_x(ig)*(N_ig'*N_ig)*dvolu;
    Kve = Kve - visc*(nx'*nx)*dvolu;
end
