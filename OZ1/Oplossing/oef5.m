e = 1e-12;
A = [1 1; 0 e; 0 0];
[Q,R] = qr(A)
b = rand(3,1);

x1 = (A'*A) \ (A'*b)
x2 = R \ (Q'*b)

norm(A*x1 - b)
norm(A*x2 - b)