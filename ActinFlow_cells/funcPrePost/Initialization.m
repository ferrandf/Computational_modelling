%CELL PARAMETERS

cellinfo = cell(1,1);
cellinfo{1}.meshparam = meshinfo(parameter);
cellinfo{1}.V_cell = zeros(parameter.nSteps,1);
tsim = 0;

%PARAMETERS FOR CONTROL PLOT AND VTK
timedata  = zeros(parameter.nSteps,1);
v_cell_info = zeros(1,parameter.nSteps);
cell_length_info = zeros(1,parameter.nSteps);
velocitydata = zeros(cellinfo{1}.meshparam.Nnodes,parameter.nSteps);
wdata = zeros(cellinfo{1}.meshparam.Nnodes,parameter.nSteps);
vpol = zeros(cellinfo{1}.meshparam.Nnodes,parameter.nSteps);
tension_actin = zeros(cellinfo{1}.meshparam.Nnodes,parameter.nSteps);

cellDensdata(1:(cellinfo{1}.meshparam.Nnodes), parameter.nSteps) = cellinfo{1}.meshparam.DOFrho_vec;
cellDensdataM(1:(cellinfo{1}.meshparam.Nnodes), parameter.nSteps) = cellinfo{1}.meshparam.DOFrho_vecM;
