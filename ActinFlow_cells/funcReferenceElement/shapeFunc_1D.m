function [N,Nx,Nxx] = shapeFunc_1D(p,z) 
% [N,Nx,Nxx] = shapeFunc_1D(p,z) 
% 1D shape Functions for the reference element [-1,1]
%    p: degree
%    z: points at which the shape functions are evaluated (column vector)
if p == 1
    N = [0.5*(1-z)   0.5*(1+z)]; 
    Nx = [-0.5*ones(size(z))   0.5*ones(size(z))]; 
    Nxx = [zeros(size(z))   zeros(size(z))]; 
elseif p == 2
    N = [0.5*z.*(z-1)   0.5*z.*(z+1)   (1+z).*(1-z)];
    Nx = [z-0.5    z+0.5    -2*z]; 
    Nxx = [ones(size(z))   ones(size(z))   -2*ones(size(z))]; 
else
    error('not considered case')
end