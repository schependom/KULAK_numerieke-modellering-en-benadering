function G = givens(A, i, j, k)
    x = A(i,k);
    y = A(j,k);
    n = sqrt(x*x + y*y);
    c = x/n;
    s = -y/n;
    G = eye(size(A,1));
    G(i,i) = c;
    G(i,j) = s;
    G(j,i) = -s;
    G(j,j) = c;
end
