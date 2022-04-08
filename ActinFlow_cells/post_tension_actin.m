function [f,fN1,fN2] = post_tension_actin(parameter,cellinfo)

%RHS: BC's + chi*int[disc_rho*dxW]
%disc_rho = N1*rho1 + N2*rho2

%---------PARAMETERS FOR ISOPARAMETRIC INTEGRATION-------------------------
p = cellinfo{1}.meshparam.p;
ngaus = cellinfo{1}.meshparam.ngaus ;
[z, w] = quadrature_1D(ngaus);
[N,Nx, Nxx] = shapeFunc_1D(p,z);

%--------DEFINITION AND BUILDING OF GLOBAL FORCE VECTOR--------------------
nElem = cellinfo{1}.meshparam.Nele;
X = cellinfo{1}.meshparam.X;
T = cellinfo{1}.meshparam.T;
adh = zeros(nElem,ngaus);

f1 = zeros(nElem,ngaus);
f2 = zeros(nElem,ngaus);

for ielem = 1 : nElem
    Te = T(ielem,:);
    Xe = X(Te,:);
    [fe1,fe2,adh(ielem,:)] = tension(Xe, Te, parameter,cellinfo,N,Nx,w);
    f1(ielem,:) = (N*fe1)';
    f2(ielem,:) = (N*fe2)';
end

fN1 = zeros(nElem+1,1);
fN2 = zeros(nElem+1,1);
f = zeros(nElem+1,1);
counter = zeros(nElem+1,1);

for i=1:nElem+1         %loop over all nodes
    for j=1:nElem       %loop over all elements
        for k=1:ngaus   %loop over 4 nodes of the element
            if T(j,k) == i
                fN1(i) = fN1(i,1) + f1(j,k);
                fN2(i) = fN2(i,1) + f2(j,k);
                counter(i) = counter(i)+1;
            end
        end
    end
    fN1(i) = fN1(i) / counter(i);
    fN2(i) = fN2(i) / counter(i);
    f(i) = fN1(i) + fN2(i);
end
end

%--------DEFINITION OF ELEMENT FORCE VECTOR------------------------------
function  [fe1,fe2,adh_f] = tension(Xe,Te, parameter,cellinfo,N,Nx,w)

nen = cellinfo{1}.meshparam.nen;
ngaus = cellinfo{1}.meshparam.ngaus;
fe1  = zeros(nen,1);
fe2  = zeros(nen,1);

chi = parameter.chi;
visc = parameter.visc;

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
    Jacob = Nxi_ig * Xe;
    dvolu = w(ig) * det(Jacob);
    res = Jacob \ Nxi_ig;
    nx = res(1,:);
    adh_f = zeros(ngaus,1);
    
    Xe_disc  = N_ig * Xe;
    rho_disc = N_ig * rho_vec';
    v_disc   = nx * v_vec';
    rho_discM = N_ig * rho_vecM';
    
    fe1(ig) = visc * v_disc;
    fe2(ig) = chi * rho_discM * rho_disc;
    
end
end