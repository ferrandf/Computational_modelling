function uNodes=computeNodesMeanValue(u,T)

nOfNodes=max(max(T));
uNodes=zeros(nOfNodes,1); counter=zeros(nOfNodes,1);
for e=1:size(T,1) %loop in elements
    Te=T(e,:);
    uNodes(Te)=uNodes(Te)+u(e);
    counter(Te)=counter(Te)+1;
end
uNodes=uNodes./counter;
