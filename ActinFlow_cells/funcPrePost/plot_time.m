
% Control plot variables

% Vcell
v_cell_info(1,parameter.step) = cellinfo{1}.V_cell(parameter.step);    
% Length of cell
cell_length_info(1,parameter.step) = abs(-cellinfo{1}.meshparam.X(1) + cellinfo{1}.meshparam.X(end)); 
% Tension
tension_actin(:,parameter.step) = f;
% Data for X cordinates of nodes
cellxcoordsdata(1:(cellinfo{1}.meshparam.Nnodes), parameter.step) = cellinfo{1}.meshparam.X(1:end);
% Rho  
cellDensdata(1:(cellinfo{1}.meshparam.Nnodes), parameter.step) = cellinfo{1}.meshparam.DOFrho_vec;
%Adhesion
cellDensdataAdh(1:(cellinfo{1}.meshparam.Nnodes),parameter.step) = cellinfo{1}.meshparam.DOFadh_vec;

%zeta:
%cellzeta(1:(cellinfo{1}.meshparam.Nnodes),parameter.step) = cellinfo{1}.meshparam.DOFzeta_vec;

% Data for cell velocity in lab frame
velocitydata(1:(cellinfo{1}.meshparam.Nnodes), parameter.step) = cellinfo{1}.meshparam.DOFv_vec(1:end);
% Data for cell velocity in cell frame
wdata(1:(cellinfo{1}.meshparam.Nnodes), parameter.step) = cellinfo{1}.meshparam.DOFw_vec(1:end);
% Polymerisation velocity
vpol(1, parameter.step) = cellinfo{1}.meshparam.vfront;
vpol(cellinfo{1}.meshparam.Nnodes, parameter.step) = cellinfo{1}.meshparam.vrear;
% RhoM
cellDensdataM(1:(cellinfo{1}.meshparam.Nnodes), parameter.step) = cellinfo{1}.meshparam.DOFrho_vecM;
% Membrane tension
memTension(1,parameter.step) = cellinfo{1}.sigma_el(parameter.step);

%---------------------------------------------------------------------
timedata(parameter.step)   = tsim;

h = abs(cellinfo{1,1}.meshparam.X(1)-cellinfo{1,1}.meshparam.X(2));
Pe = (abs(parameter.Velect) * h) / (2 * parameter.Dm);

%% plot for no water-enhacement migration 
% Update plots
if (rem(parameter.step,parameter.nplotSteps) == 0 || parameter.step == 1) && parameter.intermplotflag == 1
        figure(1); hold on;
        y = cellinfo{1}.meshparam.X';
        x = ones(size(y))*tsim;
        z = zeros(size(x));

%------COLORMAP FOR RHO----------------------------------------------------
        col_rho = cellinfo{1}.meshparam.DOFrho_vec'; 
        subplot(2,3,1);
        title('Cell movement and Actin concentration')
        xlabel('t [min]'); ylabel('X')
        surface([x;x]/60,[y;y],[z;z],[col_rho;col_rho],'facecol','no',...
            'edgecol','interp','linew',5); colorbar;
        
        col_rhoM = cellinfo{1}.meshparam.DOFrho_vecM';
        subplot(2,3,2);
        title('Cell movement and Myosin concentration')
        xlabel('t [min]'); ylabel('X')
        surface([x;x]/60,[y;y],[z;z],[col_rhoM;col_rhoM],'facecol','no',...
            'edgecol','interp','linew',5); colorbar;

%------COLORMAP FOR W----------------------------------------------------
        col_w = cellinfo{1}.meshparam.DOFv_vec';
        subplot(2,3,3);
        surface([x;x]/60,[y;y],[z;z],[col_w;col_w],'facecol','no',...
            'edgecol','interp','linew',5); colorbar;
        title('Cell movement and velocity')
        xlabel('t [min]'); ylabel('X')


        if (parameter.step > parameter.nplotSteps+1)
           subplot(2,3,4);
            plot([timedata(parameter.step-parameter.nplotSteps) timedata(parameter.step)]/60,...
                [cell_length_info(1,parameter.step-parameter.nplotSteps) cell_length_info(1,parameter.step)],'k:','linewidth',2); hold on;
            xlabel('t [min]'); ylabel('Cell length')
            
            subplot(2,3,5);
            plot([timedata(parameter.step-parameter.nplotSteps) timedata(parameter.step)]/60,...
               [v_cell_info(1,parameter.step-parameter.nplotSteps) v_cell_info(1,parameter.step)],'k:','linewidth',2); hold on;
           xlabel('t [min]'); ylabel('ALE velocity')
    
            subplot(2,3,6)
             plot([timedata(parameter.step-parameter.nplotSteps) timedata(parameter.step)]/60,...
                [memTension(1,parameter.step-parameter.nplotSteps) memTension(1,parameter.step)],'k:','linewidth',2); hold on;
            xlabel('t [min]'); ylabel('Membrane tension')
        end


    drawnow
    
  
end
    
