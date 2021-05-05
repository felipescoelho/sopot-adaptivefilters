function y = reapprox_mat(x, maxSumsfilt)

%   reapprox_max.m
%       Function to reapproximate a matrix in the spt format.
%
%
%   Syntax:
%       y = reapprox_mat(x)
%
%
%   Input Arguments:
%       . x     : Matrix to be reapproximated.
%                   3-D ARRAY (nCoefficients, nCoefficients, bitplaneDepth)
%
%
%   Output Arguments:
%       . y     : Reapproximated matrix.
%                   3-D ARRAY (nCoefficients, nCoefficients, bitplaneDepth)
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%

[nCoefficients, ~, bitplaneDepth] = size(x);
y = zeros(nCoefficients, nCoefficients, bitplaneDepth);

for i = 1:nCoefficients
    x_aux(:, :) = x(i, :, :);
    for j = 1:nCoefficients
        y(i, j, :) = reapprox_const(x_aux(j, :), maxSumsfilt);
    end
end



% EoF

