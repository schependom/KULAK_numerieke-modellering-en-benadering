A = rand(5,5);
A(2:4,1) = 0;
A(3,2) = 0;
A(4,2) = 0;
A(4,3) = 0;

spy(A)
pause

G = givens(A, 1, 5, 1);
A = G'*A;
spy(abs(A) > 1e-12);
pause

G = givens(A, 2, 5, 2);
A = G'*A;
spy(abs(A) > 1e-12);
pause

G = givens(A, 3, 5, 3);
A = G'*A;
spy(abs(A) > 1e-12);
pause

G = givens(A, 4, 5, 4);
A = G'*A;
spy(abs(A) > 1e-12);
