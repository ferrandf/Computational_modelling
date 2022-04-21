function [cellinfo] = adhesion_dynamics(cellinfo,parameter) 

%%%%% Source term constants 
Sa = parameter.Sa;
Da = parameter.Da;
gamma = parameter.gamma;

source_rhs = Sa;


[M,K,r] = transportMatrices_Adhesion(cellinfo);

[A_aux,B_aux] = timeIntegrationMatrices_adhesion(parameter,cellinfo,M,K,Da,gamma);

theta = parameter.theta;
ind = cellinfo{1}.meshparam.ind;

% Convection matrices (time-dependent)
C1 = transportMatrices_C_adhesion(cellinfo);
 
A = A_aux;
A(ind,ind) = A(ind,ind) + theta*C1 + theta*M*gamma;
B = B_aux;
B(ind,ind) = B(ind,ind) - (1-theta)*C1 - (1-theta)*M*gamma;

rhs = B*cellinfo{1}.meshparam.DOFadh_vec + r*source_rhs;

cellinfo{1}.meshparam.DOFadh_vec = A\rhs;