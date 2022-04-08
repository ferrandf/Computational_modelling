function [cellinfo] = Transport_M_single(cellinfo,parameter) 

%%%%% Source term constants 
kd = parameter.kdM;
kp = parameter.kpM;

tau_s = 1/kd;
rho_zero = kp/kd;
source_rhs = rho_zero/tau_s;

[M,K,r] = transportMatrices_MK(cellinfo);

[A_aux,B_aux] = timeIntegrationMatrices_M(parameter,cellinfo,M,K,kd,kp,parameter.diffM);

theta = parameter.theta;
ind = cellinfo{1}.meshparam.ind;

% Convection matrices (time-dependent)
C1 = transportMatrices_C_M(cellinfo);

A = A_aux;
A(ind,ind) = A(ind,ind) + theta*C1 + theta*M/tau_s*0;
B = B_aux;
B(ind,ind) = B(ind,ind) + (theta-1)*C1 + (theta-1)*M/tau_s*0;

rhs = B*cellinfo{1}.meshparam.DOFrho_vecM  + r*source_rhs*0;

cellinfo{1}.meshparam.DOFrho_vecM = A\rhs;

