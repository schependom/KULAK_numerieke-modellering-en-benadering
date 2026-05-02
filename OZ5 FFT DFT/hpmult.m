function cv = hpmult(av, bv)
% HPMULT 
%   computes the product of two numbers in high precision 
%   by using the fft
%
% INPUTS
%       av  Vector with the first number to be multiplied (flipped left-right)
%       bv  Vector with the second number to be multiplied (flipped left-right)
%     
% OUTPUTS
%       cv  Vector with the product of the two numbers (flipped left-right)
    
    % N-1 = lengte resultaat = l
    l = max(length(av), length(bv))*2;

    % Pad both vectors
    av = [av zeros(1, l - length(av))];
    bv = [bv zeros(1, l - length(bv))];

    % Compute FFT coefficients of polynomial coefficients (p100)
    Xa = fft(av);
    Xb = fft(bv);
    
    % Element-wise product of fft of both vectors
    C = Xa.*Xb;

    % Compute inverse fft of product
    cv = round(real(ifft(C)));  
    
    % Carry over digits that exceed B=10
    for i = 1:l-1
        m = mod(cv(i),10);                      % Compute remainder of ith digit after division by 10
        cv(i+1) = cv(i+1) + floor(cv(i)/10);    % Carry excess over to next digit
        cv(i) = m;                              % Replace ith digit by its value mod 10
    end
    
    % Trim zeros from result
    cv = cv(1:find(cv ~= 0, 1, 'last'));

end