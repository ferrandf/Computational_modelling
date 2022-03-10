function u=analytical_island(X)
tol = 1.e-10;
u = zeros(size(X(:,1))); %X(:,1);
u(find(abs(X(:,1)+ 0.807396892118366) < tol)) = 6;

