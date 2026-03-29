function opgave5(sigma, U, V, type)
%OPGAVE 5 Een functie die niks returnt maar gebruik maakt van opgave4
%
% Inputs
%   sigma   n diagonaalelementen van S
%   U       (m x m) orthogonale matrix
%   V       (n x n) orthogonale matrix
%   type    bepaalt c en dc

m = size(U,1); % 10
n = size(V,1);  % 7

% veronderstel dat de sigma geordend zijn van groot naar klein
I = eye(m);

% Grootste singuliere waarde \sigma_1
e1 = I(:,1); % eerste kolom: [1;0;0;0;...;0;0]

% Kleinste singuliere waarde \sigma_n
e7 = I(:,n); % tweede kolom: [0;...;0;1;0;0;0]

if type == 1 
  c = rand(10,1); % random kolom
  dc = rand(10,1);% random kolom
elseif type == 2
  c = rand * e1;  % random element in kolom
  dc = rand * e7; % random element in kolom
elseif type == 3  
  c = rand * e7;  % random element in andere kolom
  dc = rand * e1; % random element in andere kolom
elseif type == 4  
  c = rand(10,1); % random kolom
  dc = [zeros(7,1); rand(3,1)]; % 3 random laatste elementen
elseif type == 5  
  c = rand(10,1); % random kolom
  dc = [rand(7,1); zeros(3,1)]; % 7 random eerste elementen
end

% oorspronkelijk stelsel
[A, b1, x1, kapA, kapbx, eta] = opgave4(sigma, c, U, V);

% geperturbeerd stelsel via opgave4
tol = 1e-8;
[A, b2, x2] = opgave4(sigma, c + tol * dc, U, V);

% geperturbeerd stelsel via \ operator
% b2 = b1 + tol * U * dc;
% x2 = A \ b2;

rx = norm(x2-x1) / norm(x1);
rb = norm(b2-b1) / norm(b1);
kexp = rx / rb;

fprintf(1, '%.1e  %3d  %.1e %.1e %.1e %.1e %.1e \n', ...
	kapA, type, eta, kapbx, rx, rb, kexp);


