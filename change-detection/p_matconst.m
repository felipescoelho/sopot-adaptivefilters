function z = p_matconst(x, y, maxSumsfilt, maxPot)

%   p_matconst.m
%       Implements the rpoduct of a matrix and a scalar for coefficients
%       represented in spt. This product uses the function inner_p.m and is
%       necessary for the implementation of the spt_RLS filter.
%
%
%   Syntax:
%       z = p_matconst(x, y)
%
%
%   Input Arguments:
%       . x     : array representing a matrix.
%                                           3-D ARRAY(N, N, bitplaneDepth)
%       . y     : array representing a scalar.
%                                           ROW(1, bitplaneDepth)
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

[N, ~, bitplaneDepth] = size(x);

%   Computing Operation:

z = zeros(N, N, bitplaneDepth);
for i = 1:N
    x_aux(:, :) = x(i, :, :);  % getting each row at a time
    for j = 1:N
        z(i, j, :) = inner_p(x_aux(j, :), y, maxSumsfilt, maxPot);
    end
end


% EoF



