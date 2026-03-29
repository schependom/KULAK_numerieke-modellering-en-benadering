function a = afgeleide(x,f,M)
%
% AFGELEIDE 
% 
%   computes the residue of the approximation of the differential of 
%   the function with function values f evaluated in the points x when
%   approximated by the differential of the discrete least-squares polynomial
%   approximation of degree M.
%
% INPUTS
%   x     Vector holding the points in which the function values are given
%   f     Function values in the points x
%   M     Degree of the polynomial approximation
%
% OUTPUTS
%   a   Approximate differential in the points x

    % Set weights    
    w = ones(size(x)); 
    
    % Coefficients of function approximation
    c = kkb1(x,f,w,M); 

    % Coefficients of differential approximation
    cdiff = c.*(0:M)'; 
    
    % Evaluate differential and compute maximum of residue
    a = polyval(cdiff(end:-1:2),x); % Evaluate differential
    
end
