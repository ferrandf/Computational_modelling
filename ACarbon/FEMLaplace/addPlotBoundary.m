function addPlotBoundary(X,Tb)

hold on
for i=1:size(Tb,1)
    plot(X(Tb(i,:),1),X(Tb(i,:),2),'k-')
end
hold off