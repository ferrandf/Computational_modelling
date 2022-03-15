function PlotBoundary(Tboundary,X)

hold on
for i=1:length(Tboundary)
    plot(X(Tboundary(i,:),1),X(Tboundary(i,:),2),'b-','LineWidth',2)
end
hold off