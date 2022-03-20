function PlotNodalField(u,X,T)

[nelem,nen]=size(T);

switch nen
    case 3 %Linear triangles
        p=X';
        t=[T';ones(1,nelem)];
    case 4 %Linear Quadrilaterals (split in four triangles)
        k1=T(:,1); k2=T(:,2); k3=T(:,3); k4=T(:,4); k5=size(X,1)+[1:nelem]';
        XnodesInteriors=(X(k1,:)+X(k2,:)+X(k3,:)+X(k4,:))/4;
        uNodesInteriors=(u(k1)+u(k2)+u(k3)+u(k4))/4;
        u = [u;uNodesInteriors];
        p=[X;XnodesInteriors]';
        t = [k1 k2 k5; k2 k3 k5; k3 k4 k5; k4 k1 k5]';
        t=[t ; ones(1,4*nelem)];
    case 6 %Triangles p=2
        p=X';
        k1=T(:,1); k2=T(:,2); k3=T(:,3); k4=T(:,4); k5=T(:,5); k6=T(:,6);
        t =[k1 k4 k6; k4 k5 k6; k4 k2 k5; k6 k5 k3]';
        t=[t ; ones(1,4*nelem)];
    case 9 %Quadrilaterals p=2 (split in 8 triangles)
        p=X';
        k1=T(:,1); k2=T(:,2); k3=T(:,3); k4=T(:,4); k5=T(:,5); 
        k6=T(:,6); k7=T(:,7); k8=T(:,8); k9=T(:,9); 
        t =[k1 k5 k9;k1 k9 k8; k5 k2 k6; k5 k6 k9;...
            k8 k9 k7; k8 k7 k4; k9 k6 k3; k9 k3 k7]';
        t=[t ; ones(1,8*nelem)];
    otherwise
        disp('Error element not implemented')
end
pdesurf(p,t,u), colormap default

