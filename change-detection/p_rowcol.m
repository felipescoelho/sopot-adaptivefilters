function z = p_rowcol(x, y, maxSumsfilt, maxPot)

%   p_rowcol.m
%       Implements the product of a row matrix with a column matrix with
%       cofficients represented in spt. This product uses the function
%       inner_p.m and it is necessary for the implementation of the spt-RLS
%       filter.
%
%
%   Syntax:
%       z = p_rowcol(x, y)
%
%
%   Input Arguments:
%       . x     : array representing the row vector.
%                                               2-D ARRAY(bitplaneDepth, N)
%       . y     : array representing the column vector.
%                                               2-D ARRAY(N, bitplaneDepth)
%
%
%   Output Arguments:
%       . z     : array containing the resulting number.
%                                               ROW(1, bitplaneDepth)
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%

[bitplaneDepth, N] = size(x);

%   Computing Operation:

z_aux1 = zeros(N, bitplaneDepth);
for i = 1:N
    z_aux1(i, :) = inner_p(transpose(x(:,i)), y(i,:), maxSumsfilt, maxPot);
end


%   Reapproximation Process:

z_aux2 = sum(z_aux1, 1);
clear z_aux1
count = 0;

while norm(z_aux2, inf) > 1
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
       % fprintf('P_ROWCOL: overflow.\n')
    elseif abs(z_aux2(1)) == 2 && count == 0
        z_aux = sign(z_aux2(1))*ones(1, bitplaneDepth);
        z_aux2(1) = z_aux2(1) - sign(z_aux2(1));
        z_aux2 = z_aux2 + z_aux;
        count = 1;
    end
end

%   adjusting number of non-zero terms
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
z = z_aux2;


% EoF



