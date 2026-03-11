% Klassieke Gram-Schmidt
% (economy size)

m = 5;
n = 4;

A = randi(10,m,n) % 4x5 matrix


%% Met twee for-lussen

Q = zeros(m,n);
R = zeros(n);

for j = 1:n

    v_j = A(:,j);

    for i = 1:(j-1)
        q_i = Q(:,i);
        a_j = A(:,j);
        r_ij = q_i'*a_j;
        v_j = v_j - r_ij*q_i;
        R(i,j) = r_ij;
    end

    r_jj = norm(v_j, 2);
    R(j,j) = r_jj;
    Q(:,j) = v_j / r_jj;

end

Q
R

% Uniek OP HET TEKEN NA
[Qm, Rm] = qr(A, "econ")
abs(Q)-abs(Qm)
abs(R)-abs(Rm)


%% Met een for-lus

Q = zeros(m,n);
R = zeros(n);

for j = 1:n

    v_j = A(:,j);

    % We projecteren a_j (v_j) op alle vorige kolommen
    % van Q die we al geconstrueerd hebben, dus m.a.w.
    % t.e.m. kolom j-1
    proj = Q(:, 1:j-1)'*v_j;

    % Dit zijn voor elke j j-1 coefficienten die
    % moeten worden afgetrokken van de huidige v_j
    % en komen terecht in de R matrix in kolom j
    % boven de diagonaal
    R(1:j-1,j) = proj;

    % Nu moeten we op basis van deze scalars de q
    % vectoren aftrekken van a
    v_j = v_j - Q(:, 1:j-1)*proj;

    r_jj = norm(v_j, 2);
    R(j,j) = r_jj;
    Q(:,j) = v_j / r_jj;

end

[Qm, Rm] = qr(A, "econ")
abs(Q)-abs(Qm)
abs(R)-abs(Rm)