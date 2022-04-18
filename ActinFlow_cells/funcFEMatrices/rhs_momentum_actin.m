function [f,adhe] = rhs_momentum_actin(parameter,cellinfo)

%RHS: BC's + chi*int[disc_rho*dxW]
%disc_rho = N1*rho1 + N2*rho2

%---------PARAMETERS FOR ISOPARAMETRIC INTEGRATION-------------------------
p = cellinfo{1}.meshparam.p;                                        
ngaus = cellinfo{1}.meshparam.ngaus ;  
[z, w] = quadrature_1D(ngaus);
[N,Nx, Nxx] = shapeFunc_1D(p,z);

%--------DEFINITION AND BUILDING OF GLOBAL FORCE VECTOR--------------------
f = zeros(cellinfo{1}.meshparam.Nnodes,1);
nElem = cellinfo{1}.meshparam.Nele;
X = cellinfo{1}.meshparam.X;
T = cellinfo{1}.meshparam.T;
adh = zeros(nElem,ngaus);

for ielem = 1 : nElem
    Te = T(ielem,:);
    Xe = X(Te,:);
    q = [ielem ielem+1]';
    [fe,adh(ielem,:)] = rhs_ele_actin(Xe, Te, parameter,cellinfo,N,Nx,w);
    adhe(ielem) = (adh(ielem,1) + adh(ielem,2))/2;
    f(q,1) = f(q,1) + fe;
end

%--------DEFINITION OF ELEMENT FORCE VECTOR------------------------------
function  [fe,adh_f] = rhs_ele_actin(Xe,Te, parameter,cellinfo,N,Nx,w)
nen=cellinfo{1}.meshparam.nen;
ngaus = cellinfo{1}.meshparam.ngaus ; 
fe  = zeros(nen,1);
chi = parameter.chi;  

rho_vec = [cellinfo{1}.meshparam.DOFrho_vec(Te(1,1),1)...
         cellinfo{1}.meshparam.DOFrho_vec(Te(1,2),1)];

rho_vecM = [cellinfo{1}.meshparam.DOFrho_vecM(Te(1,1),1)...
         cellinfo{1}.meshparam.DOFrho_vecM(Te(1,2),1)];
     
     
     
v_vec = [cellinfo{1}.meshparam.DOFv_vec(Te(1,1),1)...
         cellinfo{1}.meshparam.DOFv_vec(Te(1,2),1)];  

%--------------BUILDING THE ELEMENT FORCE VECTOR---------------------------
for ig = 1:ngaus
    N_ig    = N(ig,:);
    Nxi_ig  = Nx(ig,:);
    Jacob = Nxi_ig*Xe;
    dvolu = w(ig)*det(Jacob);
    res = Jacob\Nxi_ig;
    nx = res(1,:);
    adh_f=zeros(ngaus,1);
    
    rho_disc = N_ig*rho_vec';
    rho_discM = N_ig*rho_vecM';

    % Multiply S to chi
    fe = fe + chi*rho_disc*rho_discM*nx'*dvolu; 
end



    













