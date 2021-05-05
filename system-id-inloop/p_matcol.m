function z = p_matcol(x, y, maxSumsfilt, maxPot)

%   p_matcol.m
%       Implements the product of a matrix and a column vector for
%       coefficients represented in spt. This product uses the function
%       inner_p.m and is necessary for the implementation of the spt-RLS
%       filter.
%
%
%   Syntax:
%       z = p_matcol(x, y)
%
%
%   Input Arguments:
%       . x     : array representing a matrix.
%                                           3-D ARRAY(N, N, bitplaneDepth)
%       . y     : array representing a column vector.
%                                           2-D ARRAY(N, bitplaneDepth)
%
%
%   Output Arguments:
%       . z     : array containing the resulting column vector.
%                                           2-D ARRAY(N, bitplaneDepth)
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%

[N, dummy1, bitplaneDepth] = size(x);
[dummy2, dummy3] = size(y);

if bitplaneDepth ~= dummy3 || (dummy2 ~= dummy1 && dummy1 ~= N)
    fprintf('The input for p_matcol.m is not in the appropriate shape')
    return
end
clear dummy1 dummy2 dummy3


%   Computing Operation:

z_aux1 = zeros(N, bitplaneDepth);
z = zeros(N, bitplaneDepth);
for k = 1:N
    x_aux(:, :) = x(k, :, :);
    for j = 1:N
        z_aux1(j, :) = inner_p(x_aux(j, :), y(j, :), maxSumsfilt, maxPot);
    end
    %   Repproximation Process:
    z_aux2 = sum(z_aux1, 1);
    count = 0;
    while norm(z_aux2, inf) > 1
        for i = 2:bitplaneDepth
            if z_aux2(i-1) > 0.0 && z_aux2(i) < 0.0
                z_aux2(i-1) = z_aux2(i-1) - 1;
                z_aux2(i) = z_aux2(i) + 2;
            elseif z_aux2(i-1) < 0.0 && z_aux2(i) > 0.0
                z_aux2(i-1) = z_aux2(i-1) + 1;
                z_aux2(i) = z_aux2(i) - 2;
            end
        end
        % adjusting for absolute values greater than 1
        for i = bitplaneDepth:-1:2
            z_aux2(i-1) = z_aux2(i-1) + sign(z_aux2(i))*floor(abs(z_aux2(i))/2);
            z_aux2(i) = z_aux2(i) - 2*sign(z_aux2(i))*floor(abs(z_aux2(i))/2);
        end
        % adjusting near spt terms
        for i = 2:bitplaneDepth
            if z_aux2(i-1) > 0.0 && z_aux2(i) < 0.0
                z_aux2(i-1) = z_aux2(i-1) - 1;
                z_aux2(i) = z_aux2(i) + 2;
            elseif z_aux2(i-1) < 0.0 && z_aux2(i) > 0.0
                z_aux2(i-1) = z_aux2(i-1) + 1;
                z_aux2(i) = z_aux2(i) - 2;
            end
        end
        % overflow
        if (abs(z_aux2(1)) >= 2 && count == 1) || (abs(z_aux2(1)) > 2)
            z_aux2 = sign(z_aux2(1))*ones(1, bitplaneDepth);
            fprintf('P_MATCOL: overflow.\n')
        elseif abs(z_aux2(1)) == 2 && count == 0
            z_aux = sign(z_aux2(1))*ones(1, bitplaneDepth);
            z_aux2(1) = z_aux2(1) - sign(z_aux2(1));
            z_aux2 = z_aux2 + z_aux;
            count = 1;
        end
    end
    %   adjusting the number of non-zero terms
        nz = nnz(z_aux2);
    while nz > maxSumsfilt
        for it = bitplaneDepth:-1:1
            if z_aux2(it) ~= 0
                z_aux2(it) = 0;
                nz = nnz(z_aux2);
                break
            end
        end
    end
    z(k, :) = z_aux2;
end

% EoF



