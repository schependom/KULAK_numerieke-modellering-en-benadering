function c = kkb2simple(x, f, w, M)
%KKB2v2 
% 
% Discrete weighted least-squares approximation of degree n for the
% function that takes the values f in the points x for the given weigths w.
%
% This method uses the normal equations to solve the system of equations.
%
% ! Solving Gc=F is equivalent to solving an overdetermined system !
%       -> see kkb2.m
%
% Inputs
%   x     Vector holding the points in which function values are given
%   f     Function values in the given points
%   w     Vector with weights in given points
%   M     Degree of the approximation
%
% Outputs
%   c     Coefficients of the approximating degree n polynomial
   
    % Make sure that all vectors are column vectors
    x = x(:); f = f(:); w = w(:);
    
    % Aantal datapunten N moet minstens gelijk zijn aan het aantal 
    % coefficienten (M+1) om een unieke oplossing te hebben.
    %
    % We laten overgedetermineerd & interpolatie toe
    %
    % N > M+1  -> Overgedetermineerd (Kleinste Kwadraten)
    % N = M+1  -> Interpolatie (Uniek bepaald, residu = 0)
    % N < M+1  -> Ondergedetermineerd (Oneindig veel oplossingen) -> ERROR
    N = length(x);
    if N < M + 1
        error('plotres:input', ...
         'length(x) should be larger than n for an overdetermined system');  
    end   
    if any([length(f),length(w)]~= length(x))
        error('plotres:input', ...
                'r and w should have the same number of elements as x');  
    end
    
    % Form normal equation explicitly and solve the system
    W = diag(w);
    A = cumprod([ones(N, 1) x(:,ones(1,M))], 2);
    b = A'*W*f; % of A'*(w.*f)
    c = (A'*W*A) \ b; % of A'*(w.*A)

end
