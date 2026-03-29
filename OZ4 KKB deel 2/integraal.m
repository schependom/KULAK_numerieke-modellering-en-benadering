function I = integraal(x,f,n)
% INTEGRAAL 
%       
%       Computes the integral of the function with function values f 
%       evaluated in the points x by computing the integral of the discrete 
%       least-squares polynomial approximation of degree n.
%
% INPUTS
%   x     Vector holding the points in which the function values are given
%   f     Function values in the points x
%   n     Degree of the polynomial approximation
%
% OUTPUTS
%   I     Value of the integral

    % Set weights
    w = ones(size(x)); 

    % Coefficients of function approximation (using QR-fact version of overdetermined system)
    c = kkb1(x,f,w,n);

    % Coefficients of integral approximation
    %   -> coefficienten van de primitieve veelterm van de gegeven veelterm
    cint = [0;c./(1:n+1)']; 
    
    % Compute integral by subtracting the values in the upper and lower limits
    Ilimits = polyval(cint(end:-1:1),[x(end),x(1)]); % Upper and lower limit
    I = Ilimits(1) - Ilimits(2); % Subtract
    
    x = linspace(0,1);
    % figure;
    cint
    % plot(x,polyval(cint(end:-1:1),x));
    
end

