function [cellinfo] = adhesion_dynamics(cellinfo,parameter) 

%%%%% Source term constants 
kd = parameter.kd;
kp = parameter.kp;

tau_s = 1/kd;
rho_zero = kp/kd;
source_rhs = rho_zero/tau_s;

[M,K,r] = transportMatrices_MK(cellinfo);

[A_aux,B_aux] = timeIntegrationMatrices_adhesion(parameter,cellinfo,M,K,kd,kp,parameter.diff);

theta = parameter.theta;
ind = cellinfo{1}.meshparam.ind;

% Convection matrices (time-dependent)
C1 = transportMatrices_C_actin(cellinfo);
 
A = A_aux;
A(ind,ind) = A(ind,ind) + theta*C1 + theta*M*gamma;
B = B_aux;
B(ind,ind) = B(ind,ind) - (1-theta)*C1 + (theta-1)*M*gamma  ;

rhs = B*cellinfo{1}.meshparam.DOFrho_vec + r*source_rhs;

cellinfo{1}.meshparam.DOFrho_vec = A\rhs;