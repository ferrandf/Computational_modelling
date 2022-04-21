clear; close all; clc
tic

%% Paths
restoredefaultpath
addpath('funcReferenceElement')
addpath('funcFEMatrices')
addpath('funcPrePost')

%% Set flags and parameters for the simulation
Parameters;
Initialization;

%% Main loop
for n = 1:parameter.nSteps
    parameter.step = n;
    
    %Momentum balance
    [cellinfo] = momentum_balance_actin(parameter,cellinfo); %solve deltax sigma = mu v (v = velocity in lab frame is plotted in figure 1 chart 3)
    
    %Ale update
    [cellinfo] = ALE_meshvel(cellinfo, parameter); %change in the coordinate system. (compute the velocity of the cell and w = vel at the cell frame) because
    % in the convective term, we need w, not v. (w = v-vcell)
    
    %Model
    [cellinfo] = Transport_actin_single(cellinfo,parameter); %computes density of actin: rhoF
    [cellinfo] = Transport_M_single(cellinfo,parameter); %computes density of myosin: rhoM. together we have K rho = f 
    [cellinfo] = adhesion_dynamics(cellinfo,parameter); %adhesion term.
    %[cellinfo] = zeta(cellinfo,parameter); %zeta
    
    %Tension
    [f,~,~] = post_tension_actin(parameter,cellinfo);

    %Update simulation variables
    tsim = tsim + parameter.Dt;
    cellinfo{1}.meshparam.X = cellinfo{1}.meshparam.X + ...
        cellinfo{1}.meshparam.DOFv_ale_vec * parameter.Dt;
    
    %On-line plots
    plot_time;
end

%% Post-process
Control_plot;

toc