function cellinfo = ALE_meshvel(cellinfo,parameter)

%  compute v front
% V0=0.55um/s from Mogilner stimates
% V=0.66-0.33 um/s from Guathier stimates
% dekta/kBT = 4.10-3/4.1pn/Um = 1kPa

n = parameter.step;
V0 = parameter.V0;
coeff1 = parameter.coeff1;

coeff3 = 4;

if cellinfo{1}.sigma_el(parameter.step)<coeff1
    cellinfo{1}.meshparam.vfront =  (V0*(1-cellinfo{1}.sigma_el(parameter.step)/coeff1)^coeff3);
    cellinfo{1}.meshparam.vrear  = (-V0*(1-cellinfo{1}.sigma_el(parameter.step)/coeff1)^coeff3);
else 
    cellinfo{1}.meshparam.vfront = 0;
    cellinfo{1}.meshparam.vrear  = 0; 
end        

% Update ALE mesh velocity - Only one density
cellinfo{1}.meshparam.DOFv_ale_vec(1) =  ...
    cellinfo{1}.meshparam.DOFv_vec(1)+cellinfo{1}.meshparam.vrear;
cellinfo{1}.meshparam.DOFv_ale_vec(end) = ...
    cellinfo{1}.meshparam.DOFv_vec(end)+cellinfo{1}.meshparam.vfront;

% Linear interpolation for midnodes
cellinfo{1}.meshparam.DOFv_ale_vec(2:end-1) = cellinfo{1}.meshparam.DOFv_ale_vec(1) +...
    ( cellinfo{1}.meshparam.X(2:end-1) - cellinfo{1}.meshparam.X(1))...
    *  (cellinfo{1}.meshparam.DOFv_ale_vec(end) - cellinfo{1}.meshparam.DOFv_ale_vec(1))...
    / ( cellinfo{1}.meshparam.X(end) - cellinfo{1}.meshparam.X(1));

% CALCULATING VELOCITY IN CELL FRAME w: If static cell, w and v are equal

cellinfo{1}.meshparam.DOFw_vec(1:end) = ...
    cellinfo{1}.meshparam.DOFv_vec(1:end)...
    - cellinfo{1}.meshparam.DOFv_ale_vec(1:end);

%Cell velocity Vcell
cellinfo{1}.V_cell(parameter.step) = (cellinfo{1}.meshparam.DOFv_ale_vec(1) + ...
    cellinfo{1}.meshparam.DOFv_ale_vec(end))/2.0;

end

