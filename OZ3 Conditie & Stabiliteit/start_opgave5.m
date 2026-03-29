m = 10;
n = 7;

% Orthogonale (vierkante) matrices
U = orth(rand(m));
V = orth(rand(n));

% We genereren achtereenvolgens problemen met
% kappa(A) = 1,10^3,10^6

% Header van de resultatentabel
fprintf(1, '%7s  %4s %7s %7s %7s %7s %7s \n', ...
	'kap_A', 'type', 'eta', 'kap_bx', 'rx', 'rb', 'k_exp');

% Loop over verschillende condities voor A
for kappa = [1 1e3 1e6]

  % Maak -> kolomvector (')
  %      -> n diagonaalelementen 
  %      -> tussen kappa en 1
  %      -> van DALENDE waarde
  sigma = linspace(kappa, 1, n)'; % transpose -> kolom!

  % Voor elk type, roep functie aan
  %      -> deze functie print zelf!
  for type = [1 2 3 4 5]
    opgave5(sigma, U, V, type);
  end

  % Print nieuwe lijn van tabel
  fprintf('\n')
  
end

