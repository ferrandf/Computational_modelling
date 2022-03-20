function theReferenceElement = createReferenceElemeny(degree,typeOfElement)
%theReferenceElement = createReferenceElement(degree,typeOfElement)
%Creates a struct with the information of the reference element, for
%triangle discretizations

switch degree
    case 1
        if typeOfElement==1 %TRI
            z = [0.5 0; 0 0.5; 0.5 0.5];
            w = [1 1 1]/6;
            xi = z(:,1); eta=z(:,2);
            N = [1-xi-eta,   xi,   eta];
            Nxi  = [-ones(size(xi))   ones(size(xi))    zeros(size(xi))];
            Neta = [-ones(size(xi))   zeros(size(xi))   ones(size(xi))];
            nodesCoord = [0,0;1,0;0,1];
        else
            z1=sqrt(3)/3;
            z=[-z1 -z1;-z1 z1;z1 -z1; z1 z1];
            w=[1 1 1 1];
            xi = z(:,1); eta=z(:,2);
            N= [(1-xi).*(1-eta)/4, (1+xi).*(1-eta)/4, ...
                (1+xi).*(1+eta)/4, (1-xi).*(1+eta)/4];
            Nxi= [(eta-1)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4];
            Neta= [(xi-1)/4, -(1+xi)/4,   (1+xi)/4,  (1-xi)/4 ];
            nodesCoord = [-1,-1;1,-1;1,1;-1,1];
        end
    otherwise
        error('Element not implemented yet')
end
if typeOfElement==1
    typeElement='TRI';
else
    typeElement='QUA';
end
theReferenceElement=struct('IPweights',w,'IPcoord',z,'N',N,'Nxi',Nxi,'Neta',Neta,'type',typeElement,'nodesCoord',nodesCoord);
