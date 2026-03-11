%%% Oefening 6

A = rand(100,50);

disp('Twee for-loops');
tic
[Q,R] = qr2for(A);
toc 
fprintf('Elapsed time with timeit is: %d seconds \n\n', timeit(@() qr2for(A)));
fprintf('Norm of error:\n')
disp(norm(A - Q*R));
pause;


disp('Eén for-loop');
tic
[Q,R] = qr1for(A);
toc
fprintf('Elapsed time with timeit is: %d seconds \n\n', timeit(@() qr1for(A)));
fprintf('Norm of error:\n')
disp(norm(A - Q*R));
pause;

disp('Eén for-loop, verbeterd');
tic
[Q,R] = qrtop(A);
toc 
fprintf('Elapsed time with timeit is: %d seconds \n\n', timeit(@() qrtop(A)));
fprintf('Norm of error:\n')
disp(norm(A - Q*R));
pause;

disp('Ingebouwde functie');
tic
[Q,R] = qr(A);
toc
fprintf('Elapsed time with timeit is: %d seconds \n\n', timeit(@() qr(A)));
fprintf('Norm of error:\n')
disp(norm(A - Q*R));
pause;

disp('Ingebouwde functie, gereduceerd');
tic
[Q,R] = qr(A,0);
toc
fprintf('Elapsed time with timeit is: %d seconds \n\n', timeit(@() qr(A,0)));
fprintf('Norm of error:\n')
disp(norm(A - Q*R));

%% Functies

function [Q,R] = qr2for(A)
    Q = zeros(size(A));
    R = zeros(size(A,2));
    for j = 1:size(A,2)
        v = A(:,j);
        for i = 1:j-1
            R(i,j) = Q(:,i)'*A(:,j);
            v = v - Q(:,i)*R(i,j);
        end
        R(j,j) = norm(v);
        Q(:,j) = v/R(j,j);
    end
end

function [Q,R] = qr1for(A)
    Q = zeros(size(A));
    R = zeros(size(A,2));
    
    % First column 
    v = A(:,1);
    R(1,1) = norm(v);
    Q(:,1) = v/R(1,1);
        
    % Other columns
    for j = 2:size(A,2)
        tmp = Q(:,1:j-1);
        R(1:j-1,j) = A(:,j)'*conj(tmp);
        v = A(:,j) - tmp*R(1:j-1,j);
        
        R(j,j) = norm(v);
        Q(:,j) = v/R(j,j);
    end
end

function [Q,R] = qrtop(A)
    Q = zeros(size(A));
    R = zeros(size(A,2));
    
    % First column 
    v = A(:,1);
    R(1,1) = norm(v);
    Q(:,1) = v/R(1,1);
        
    % Other columns
    for j = 2:size(A,2)
        R(:,j) = A(:,j).'*conj(Q);
        v = A(:,j) - Q*R(:,j);
        
        R(j,j) = norm(v);
        Q(:,j) = v/R(j,j);
    end
end