function z = p_colrow(x, y, maxSumsfilt, maxPot)

%   p_colrow.m
%       Implements the product of a column vector and a row vector with
%       coefficients represented in spt. This function uses the inner_p.m
%       function and is necessary for the implementation of the spr-RLS
%       filter.
%
%
%   Syntax:
%       z = p_colrow(x, y)
%
%
%   Input Arguments:
%       . x     : array representing the column vector.
%                                               2-D ARRAY(N, bitplaneDepth)
%       . y     : array representing the row vector.
%                                               2-D ARRAY(bitplaneDepth, N)
%
%
%   Output Arguments:
%       . z     : array containing the resulting matrix.
%                                           3-D ARRAY(N, N, bitplaneDepth)
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%

[N, bitplaneDepth] = size(x);
[dummy1, dummy2] = size(y);

if bitplaneDepth ~= dummy1 || N ~= dummy2
    fprintf('The input for p_colrow.m is not in the appropriate shape')
    return
end
clear dummy1 dummy2


%   Computing Operation:

z = zeros(N, N, bitplaneDepth);
for i = 1:N
    for j = 1:N
        z(i, j, :) = inner_p(x(i, :), transpose(y(:, j)), maxSumsfilt,...
            maxPot);
    end
end


% EoF



