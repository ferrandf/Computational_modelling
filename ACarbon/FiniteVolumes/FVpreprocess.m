function [elementsForSide,sidesForElement,Tsides,nOfInteriorSides,Ve,Ls,normalVectors,XmidVolumes,XmidSides] = FVpreprocess(X,T)

nOfElementNodes=size(T,1);
n=sqrt(nOfElementNodes);
if round(n)==n %quadrilatera
    n = 4;
else % triangles
    n = 3;
end
[intSides, extSides] = GetSides(T(:,1:n));

nOfElements = size(T,1);
nOfInteriorSides = size(intSides,1);
nOfExteriorSides = size(extSides,1);

F = zeros(nOfElements,3);
for iFace = 1:nOfInteriorSides
    infoFace = intSides(iFace,:);
    F(infoFace(1),infoFace(2)) = iFace;
    F(infoFace(3),infoFace(4)) = iFace;
end

for iFace = 1:nOfExteriorSides
    infoFace = extSides(iFace,:);
    F(infoFace(1),infoFace(2)) = iFace + nOfInteriorSides;
end

sidesForElement=F;
infoSides.intSides = intSides;
infoSides.extSides = extSides;
%Ghost elements are included (one per exterior side)
elementsForSide=[infoSides.intSides(:,[1,3]);...
    [infoSides.extSides(:,1),size(T,1)+[1:nOfExteriorSides]']];

[Ve,Ls,normalVectors,XmidSides,Tsides]=computeSizesAndNormals(X,T,infoSides);

XmidVolumes=(X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:))/3;

end

%______________________________________________
function [Ve,Ls,normalVectors,XmidSides,Tsides]=computeSizesAndNormals(X,T,infoSides)
% [Ve,Le,normalVectors]=computeSizesAndNormals(X,T,infoSides)
% Only for TRIANGLES

nOfVolumes=size(T,1);
Ve=zeros(1,nOfVolumes);
for e=1:nOfVolumes
    Xe=X(T(e,:),:);
    Ve(e)=det([Xe(2,:)-Xe(1,:);Xe(3,:)-Xe(1,:)])/2;
end

allSides=[infoSides.intSides(:,1:2);infoSides.extSides];
nOfSides=size(allSides,1);
nodesSides=[1 2;2 3; 3 1];
Ls=zeros(1,nOfSides); normalVectors=zeros(2,nOfSides);
XmidSides=zeros(nOfSides,2);
Tf=zeros(nOfSides,2);
for s=1:nOfSides
    Tf=T(allSides(s,1),nodesSides(allSides(s,2),:));
    Tsides(s,:)=Tf;
    Xf=X(Tf,:);
    XmidSides(s,:)=(Xf(1,:)+Xf(2,:))/2;
    t=Xf(2,:)-Xf(1,:);
    Ls(s)=norm(t);
    normalVectors(:,s)=[t(2);-t(1)]/Ls(s);
end

end
%______________________________________________
function [intSides,extSides] = GetSides(T)
%
% [intSides,extSides] = GetSides(T)
% (only for triangles and tetrahedrons)
%
% For every face i:
% intSides(i,:)=[element1 nface1 element2 nface2 node1] for interior sides
% extSides(i,:)=[element1 nface1] for exterior sides
%
% element1, element2:   number of the elements
% nface1, nface2:       number of face in each element
% node1:  number of node in the 2nd element that matches with the 1st node
%         in the 1st element (with the numbering of the face)
%
% Input:
% T: connectivity of the mesh
%
% Output:
% intSides,extSides: interior and exterior sides
%

[nElem,nen]=size(T);
nfaceel = nen;

nNodes = max(max(T));
N = zeros(nNodes,10);
nn = ones(nNodes,1);
for ielem = 1:size(T,1)
    Te = T(ielem,:);
    nn_Te = nn(Te);
    for kk = 1:nfaceel
        N(Te(kk),nn_Te(kk)) = ielem;
    end
    nn(Te) = nn(Te) +1;
end
N(:,max(nn):end) = [];

markE = zeros(nElem,nfaceel);
intSides = zeros(fix(3/2*size(T,1)),5);
extSides = zeros(size(T,1),2);

%Definition of the sides in the reference element
switch nen
    case 3 %triangle
        Esides = [1 2; 2 3; 3 1];
    case 4 %tetrahedra
        Esides = [1 2; 2 3; 3 4; 4 1];
end
intF = 1;
extF = 1;
for iElem=1:nElem
    for iFace=1:nfaceel
        if(markE(iElem,iFace)==0)
            markE(iElem,iFace)=1;
            nodesf = T(iElem,Esides(iFace,:));
            
            jelem = FindElem(iElem,nodesf,N,T);
            
            if(jelem~=0)
                [jface,node1]=FindFace(nodesf,T(jelem,:),Esides);
                intSides(intF,:)=[iElem,iFace,jelem,jface,node1];
                intF = intF +1;
                markE(jelem,jface)=1;
            else
                extSides(extF,:)=[iElem,iFace];
                extF = extF + 1;
            end
        end
    end
end

intSides = intSides(intSides(:,1)~=0,:);
extSides = extSides(extSides(:,1)~=0,:);

end

%Auxiliar functions
function jelem = FindElem(iElem,nodesf,N,T)

nen = length(nodesf);

% [elems,aux] = find(T==nodesf(1));
elems = N(nodesf(1),(N(nodesf(1),:)~=0));
elems=elems(elems~=iElem);
Ti=T(elems,:);
for i=2:nen
    if(~isempty(elems))
        [aux,aux2] = find(Ti==nodesf(i));
        elems = elems(aux);
        Ti=Ti(aux,:);
    end
end

if(isempty(elems))
    jelem=0;
else
    jelem=elems(1);
end

end

%_______________
function [jface,node1]=FindFace(nodesf,nodesE,Esides)

nSides = size(Esides,1);
for j=1:nSides
    nodesj = nodesE(Esides(j,:));
    if (nodesj(1)==nodesf(1)|| nodesj(1)==nodesf(2)) && ...
            (nodesj(2)==nodesf(1)|| nodesj(2)==nodesf(2))
        jface = j;
        node1 = find(nodesj==(nodesf(1)));
        break;
    end
end

end


