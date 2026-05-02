function n = vec2str(vecn)
% VEC2STR 
%   Convert a character array of digits, flipped left to right, to a
%   string representation of the number.
%
% INPUTS
%       vecn Character array of digits, flipped left to right
% OUTPUTS
%       n String representation of the number

    % Flip and convert character array to string
    n = num2str(fliplr(vecn)')';

end