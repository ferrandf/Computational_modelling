function [A,B] = timeIntegrationMatrices_adhesion(parameter,cellinfo,M,K,Da,gamma)

theta = parameter.theta;
Dt = parameter.Dt;
ndof = cellinfo{1}.meshparam.DOFadh;
nNodes = parameter.nx+1;
ind = cellinfo{1}.meshparam.ind;

M_dt = M/Dt; 

nNonZero = nNodes*3*5; 
A = spalloc(ndof,ndof,nNonZero); 
B = spalloc(ndof,ndof,nNonZero); 

aux_11 = Da*K ; 

A(ind,ind) = M_dt + theta*aux_11; 
B(ind,ind) = M_dt - (1-theta)*aux_11; 




