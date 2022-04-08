function [cellinfo] = momentum_balance_actin(parameter, cellinfo)

%-- -2*mu*int[dxWdxV] - eta*int[V*W] = BC's + chi*int[rho*dxW]
[K] = momentumMatrice_K_actin(parameter,cellinfo);

%-- Force vector f_bd (bd = BEFORE PRESSURE)
[f_bp,adh] = rhs_momentum_actin(parameter,cellinfo);

cellinfo{1}.meshparam.adh = adh;

f = f_bp;

L_current = abs(cellinfo{1}.meshparam.X(end,1) - ...
    cellinfo{1}.meshparam.X(1,1));
    
cellinfo{1}.sigma_el(parameter.step) = parameter.k * ...
    (L_current - parameter.L0 - parameter.Lpre)/(parameter.L0 + parameter.Lpre);


    DOFv_vec = K\f;
    cellinfo{1}.meshparam.DOFv_vec = DOFv_vec;
end