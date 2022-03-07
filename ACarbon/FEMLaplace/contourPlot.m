function contourPlot(u,X,T)

switch size(T,2)
    case 3 %linear triangles
        p=X';
        t=[T';ones(1,size(T,1))];
    case 4 %linear quadrilateral
        p=X';
        t=[[T(:,[1,2,3]);T(:,[1,3,4])]';ones(1,2*size(T,1))];       
    otherwise
        error('This function work only for linear FEM :-(')
end
  
pdecont(p,t,u,30);