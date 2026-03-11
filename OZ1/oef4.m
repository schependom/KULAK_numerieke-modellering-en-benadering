X = rand(5,2); 
X(2:5,1) = X(2:5,1)*1.e-8

u_1 = X(:,1) + norm(X(:,1), 2)*ones(5)
u_2 = X(:,1) - norm(X(:,1), 2)*ones(5)

% Twee bijna even grote getallen van elkaar aftrekken