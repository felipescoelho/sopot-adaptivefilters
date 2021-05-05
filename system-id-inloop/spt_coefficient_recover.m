function y = spt_coefficient_recover(x, maxPot)

%   spt_coefficient_recover.m
%       This function is responsable for the recovering of the filter's spt
%       approximated coefficients.
%
%
%   Syntax:
%       y = spt_coefficient_recover(x)
%
%
%   Input Arguments:
%       . x     : coefficients from the trained filter.
%           3-D ARRAY (nCoefficients, bitplaneDepth, (nIterations + 1))
%
%
%   Output Arguments:
%       . y     : recovered coefficients.
%                           2-D ARRAY (nCoefficients, (nIterations + 1))
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%

[nCoefficients, ~, nIterations_plus1] = size(x);
y = zeros(nCoefficients, nIterations_plus1);

for i = 1:nCoefficients
    x_aux1(:, :) = x(i, :, :);
    x_aux2 = x_aux1.';  % must be [nCoefficients, bitplaneDepth]
    y(i, :) = spt_signal_recover(x_aux2, maxPot);
end



% EoF

