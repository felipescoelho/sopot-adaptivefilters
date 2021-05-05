function y = reapprox_col(x, maxSumsfilt)

%   reapprox_col.m
%       The function reapproximate a spt column vector.
%
%
%   Syntax:
%       y = reapprox_col(x)
%
%
%   Input Arguments:
%       . x     : columns vector in the spt approximation format.
%                               2-D ARRAY (nCoefficients, bitplaneDepth)
%
%
%   Output Arguments:
%       . y     : column vector reapproximated to the spt format.
%                               2-D ARRAY (nCoefficients, bitplaneDepth)
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%

[nCoefficients, bitplaneDepth] = size(x);
y = zeros(nCoefficients, bitplaneDepth);

for i = 1:nCoefficients
    y(i, :) = reapprox_const(x(i, :), maxSumsfilt);
end



% EoF

