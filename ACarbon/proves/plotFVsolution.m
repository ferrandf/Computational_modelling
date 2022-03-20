function plotFVsolution(X,T,u,figureNumber)

figure(figureNumber)

clf, hold on
for e=1:size(T,1)
    Te=T(e,:);
    patch(X(Te,1),X(Te,2),u(e)*[1;1;1],u(e))
end
hold off
view(3)