function [] = plot_vtk(Cname, timedata, cellxcoordsdata, cellDensdata, cellDensdataM, RSDensdata, cellDensdataG, cellDensdatam, wdata, parameter, cellinfo)
[Tdata,Xdata] = meshgrid(timedata/60,cellxcoordsdata(:,1));
nodes=[Tdata(:) Xdata(:) ];
nnodes=size(nodes,1);

nx=size(Xdata,1);
ny=size(Xdata,2);

Nx=nx-1; Ny=ny-1; 

nelem=0;
for iy=1:Ny
    for ix=1:Nx
        nelem=nelem+1;
        elem(nelem,:)=...
            [nx*(iy-1)+ix ...
            nx*(iy-1)+ix+1 ...
            nx*(iy-1)+nx+ix+1 ...
            nx*(iy-1)+nx+ix ...
            nx*(iy-1)+ix ...
            nx*(iy-1)+ix+1 ...
            nx*(iy-1)+nx+ix+1 ...
            nx*(iy-1)+nx+ix];
    end
end
     
 	% printing heading to file
f=fopen(char(Cname+string('.vtk')),'w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'Sensitivity Analysis migration\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnodes ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

	% printing coordinates
for i1=1:nnodes
    fprintf(f,'%14.8E %14.8E %14.8E\n',[nodes(i1,1:2) 0]);
    id(i1)=nodes(i1,1);
end
fprintf(f,'\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

% printing connectivity
fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*5);
for i1=1:nelem     
    fprintf(f,'%4i  %10i  %10i  %10i %10i\n',4,...
        elem(i1,1)-1,elem(i1,2)-1,elem(i1,3)-1,elem(i1,4)-1);
end
fprintf(f,'\n');
fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
for i1=1:nelem
    fprintf(f,' %4i ', 9);
end
fprintf(f,'\n'); 
% 
fprintf(f,'%s %8i\n','POINT_DATA', nnodes);
fprintf(f,'SCALARS Actin float\n');
fprintf(f,'LOOKUP_TABLE default\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
denXXX=cellDensdata(:);
for i1=1:nnodes
        fprintf(f,'%14.8E\n', denXXX(i1,1) );
end

fprintf(f,'SCALARS Myosin float\n');
fprintf(f,'LOOKUP_TABLE default\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
denX=cellDensdataM(:);
for i1=1:nnodes
        fprintf(f,'%14.8E\n', denX(i1,1) );
end
if parameter.signaling
    fprintf(f,'SCALARS rhoR float\n');
    fprintf(f,'LOOKUP_TABLE default\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
    denRS=RSDensdata(1:cellinfo{1}.meshparam.Nnodes,:);
    denR = denRS(:);
    for i1=1:nnodes
            fprintf(f,'%14.8E\n', denR(i1,1) );
    end
    
    fprintf(f,'SCALARS rhoS float\n');
    fprintf(f,'LOOKUP_TABLE default\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
    denSR=RSDensdata(cellinfo{1}.meshparam.Nnodes+1:end,:);
    denS=denSR(:);
    for i1=1:nnodes
            fprintf(f,'%14.8E\n', denS(i1,1) );
    end
end
if parameter.coupled_actomyosin
    fprintf(f,'SCALARS ActinMonomer float\n');
    fprintf(f,'LOOKUP_TABLE default\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
    denG=cellDensdataG(:);
    for i1=1:nnodes
            fprintf(f,'%14.8E\n', denG(i1,1) );
    end

    fprintf(f,'SCALARS FreeMyosin float\n');
    fprintf(f,'LOOKUP_TABLE default\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
    denm=cellDensdatam(:);
    for i1=1:nnodes
            fprintf(f,'%14.8E\n', denm(i1,1) );
    end
end    
fprintf(f,'SCALARS FlowV float\n');
fprintf(f,'LOOKUP_TABLE default\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
wXXX=wdata(:);
for i1=1:nnodes
        fprintf(f,'%14.8E\n', wXXX(i1,1) );
end

fprintf(f,'SCALARS tension float\n');
fprintf(f,'LOOKUP_TABLE default\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
tXXX=wdata(:);
for i1=1:nnodes
        fprintf(f,'%14.8E\n', tXXX(i1,1) );
end

fclose(f);