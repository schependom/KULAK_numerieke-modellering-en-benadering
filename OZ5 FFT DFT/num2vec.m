function vecn = num2vec(n)
% NUM2VEC 
%   Convert a number or string representation of a number to a
%   character array of digits, flipped left to right.
% 
% INPUTS
%       n Number or string representation of a number to convert
% 
% OUTPUTS
%       vecn Character array of digits, flipped left to right

    % Convert number to string if necessary
    if ~ischar(n)
        n = num2str(n);
    end
    
    % Convert string to character array of digits
    vecn = str2num(n')'; % Rotate string, split, and rotate back
    vecn = fliplr(vecn); % Flip string left to right

    % VOORBEELDEN
    %
    %   123 -> "123"    via num2str
    %       -> "1"
    %          "2"      character array via ' operator
    %          "3"
    %       ->  1
    %           2       num array via str2num
    %           3 
    %       -> [1 2 3]  transpose
    %       -> [3 2 1]  fliplr

end