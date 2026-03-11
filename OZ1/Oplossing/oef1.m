
A = diag(randn(1,7));
A = A+diag(rand(1,6),1);
A = A+diag(rand(1,6),-1);
A = A+diag(rand(1,5),-2);

[Q,R] = qr(A)