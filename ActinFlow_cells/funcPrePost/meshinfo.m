function [meshout] = meshinfo(parameter)

%Computational mesh parameters

meshout.Nnodes = parameter.nx + 1;  %Nodes
meshout.nen = 2;                    %Element per node
meshout.Nele = parameter.nx;        %Number of elements
meshout.X = linspace(parameter.dom(1),parameter.dom(2),meshout.Nnodes)';
meshout.T = [1:parameter.nx; 2:parameter.nx+1]';
meshout.le = abs(meshout.X(1) - meshout.X(2));      %Length of the system            
meshout.node0 = find(abs(meshout.X - parameter.dom(1))< 1e-6);
meshout.node1 = find(abs(meshout.X - parameter.dom(2))< 1e-6);
meshout.ind = 1:meshout.Nnodes;     %Index
meshout.indR = meshout.ind;         %R-index
meshout.indS = meshout.Nnodes+1:2*meshout.Nnodes; %S-index
meshout.p = 1;                      %Polynomial degree                             
meshout.ngaus = 2 ;                 %Number of gaus points

%-------------------------------------------------------------------------
% Definition of degrees of freedom

meshout.DOFv = meshout.Nnodes;      %velocity from lab frame
meshout.DOFw = meshout.Nnodes;      %velocity from cell frame
meshout.DOFv_ale = meshout.Nnodes;  %ALE velocity
meshout.DOFrho = meshout.Nnodes;    %tranport of rho

meshout.DOFvM = meshout.Nnodes;     %velocity from lab frame
meshout.DOFwM = meshout.Nnodes;     %velocity from cell frame
meshout.DOFv_aleM = meshout.Nnodes; %ALE velocity
meshout.DOFrhoM = meshout.Nnodes;   %tranport of rho

%------------------------------------------------------------------------
% Definition of solution vectors at current time-> FOR PLOTTING
meshout.DOFv_ale_vec =  zeros(meshout.DOFv_ale,1);
meshout.DOFv_vec = zeros(meshout.DOFv,1);
meshout.DOFw_vec =  zeros(meshout.DOFw,1);

meshout.DOFrho_vec =  zeros(meshout.DOFrho,1);
meshout.DOFrho_vecM =  zeros(meshout.DOFrhoM,1);

%--------------------------------------------------------------------
% Definition of solution vectors at previous time-> FOR PLOTTING
meshout.DOFv_ale_old = meshout.DOFv_ale_vec;
meshout.DOFv_old = meshout.DOFv_vec;
meshout.DOFw_old = meshout.DOFw_vec;

meshout.DOFrho_old = meshout.DOFrho;
meshout.DOFrho_oldM = meshout.DOFrhoM;

%--------------------------------------------------------------------
% Initial condition for rho

cmean    = 1;
rhoprofile = 0.9 + 0.200 * 1 * (rand(length(meshout.X),1)/2 -1);

rhoprofile = rhoprofile/mean(rhoprofile) * cmean;
rhoprofile = rhoprofile * 0 + 1.;

meshout.DOFrho_vec(1:meshout.DOFrho) = rhoprofile;
meshout.DOFrho_vecM(1:meshout.DOFrhoM) = rhoprofile;


%--------------------------------------------------------------------
% Charge data

random_ini = 0.9 + 0.2 * rand(size(meshout.DOFrho_vec));
meshout.DOFrho_p = random_ini;




end

