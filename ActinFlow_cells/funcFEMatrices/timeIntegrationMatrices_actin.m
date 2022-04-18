function [A,B] = timeIntegrationMatrices_actin(parameter,cellinfo,M,K,kd,kp,diff)

theta = parameter.theta;
Dt = parameter.Dt;
ndof = cellinfo{1}.meshparam.DOFrho;
nNodes = parameter.nx+1;
ind = cellinfo{1}.meshparam.ind;

M_dt = M/Dt; 

nNonZero = nNodes*3*5; 
A = spalloc(ndof,ndof,nNonZero); 
B = spalloc(ndof,ndof,nNonZero); 

aux_11 = diff*K ; 

A(ind,ind) = M_dt + theta*aux_11; 
B(ind,ind) = M_dt - (1-theta)*aux_11; 




