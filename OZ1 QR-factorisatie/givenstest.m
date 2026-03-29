A = rand(5,5)
G = givens(A, 1, 5, 5);
G'*A

A = rand(5,5);
A(3,2) = 0;
A(4,2) = 0;
A(4,3) = 0;

ex = 1;

%%
if ex == 1
    A = triu(randi(10,5,4));
    A = [A(:,1:2) randi(10,5,1) A(:,3:end)]

    A1 = A
    spy(abs(A1)>1e-12)
    pause
    G1 = givens(A1, 3, 4, 3);
    A1 = G1'*A1
    spy(abs(A1)>1e-12)
    pause
    G2 = givens(A1, 3, 5, 3);
    A1 = G2'*A1
    spy(abs(A1)>1e-12)

    pause
    A2 = A
    spy(abs(A2)>1e-12)
    pause
    G1 = givens(A2, 4, 5, 3);
    A2 = G1'*A2
    spy(abs(A2)>1e-12)
    pause
    G2 = givens(A2, 3, 4, 3);
    A2 = G2'*A2
    spy(abs(A2)>1e-12)
end

%% remove column
if ex == 2
    A = triu(randi(10,5,5));
    A = [A(:,1:2) A(:,4:end)]

    A1 = A
    spy(abs(A1)>1e-12)
    pause
    G1 = givens(A1, 3, 4, 3);
    A1 = G1'*A1
    spy(abs(A1)>1e-12)
    pause
    G2 = givens(A1, 4, 5, 4);
    A1 = G2'*A1
    spy(abs(A1)>1e-12)
end

%% add row
if ex == 3
    A = rand(5,5);
    [Q, R] = qr(A);
    A2 = [ rand(1,5); A ];
    Q2 = [1 0 0 0 0 0; zeros(5,1), Q];
    R2 = Q2' * A2;

    spy(abs(R2)>1e-12)
    pause
    G = givens(R2, 1, 2, 1);
    R2 = G'*R2
    spy(abs(R2)>1e-12)
    pause
    G = givens(R2, 2, 3, 2);
    R2 = G'*R2
    spy(abs(R2)>1e-12)
    pause
    G = givens(R2, 3, 4, 3);
    R2 = G'*R2
    spy(abs(R2)>1e-12)
    pause
    G = givens(R2, 4, 5, 4);
    R2 = G'*R2
    spy(abs(R2)>1e-12)
    pause
    G = givens(R2, 5, 6, 5);
    R2 = G'*R2
    spy(abs(R2)>1e-12)
end


