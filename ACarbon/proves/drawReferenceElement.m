function drawReferenteElement(theReferenceElement)

if theReferenceElement.type == 'TRI' %triangle
    plot([0,1,0,0],[0,0,1,0],'k-','LineWidth',2,'MarkerSize',12)
else %quadrilateral
    plot([-1,1,1,-1,-1],[-1,-1,1,1,-1],'ko-')
end
pgauss = theReferenceElement.IPcoord;
nodes = theReferenceElement.nodesCoord;
hold on, plot(pgauss(:,1),pgauss(:,2),'r*','LineWidth',2,'MarkerSize',12), 
plot(nodes(:,1),nodes(:,2),'ko','LineWidth',2,'MarkerSize',12),hold off
axis equal, axis off


