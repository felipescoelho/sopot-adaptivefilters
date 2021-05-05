function y = reapprox_const(x, maxSumsfilt)

%   reapprox_const.m
%       This function reshapes a constant represented in spt indexing, in
%       order for it to fit the prestablished parameters.
%
%
%   Syntax:
%       y = reapprox_const(x);
%
%
%   Input Arguments:
%       . x     : constant represented by the spt indexing.
%                                                   ROW (1, bitplaneDepth)
%
%
%   Output Arguments:
%       . y     : constant reapproximated with the spt paremeters.
%                                                   ROW (1, bitplaneDepth)
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%

bitplaneDepth = length(x);
count = 0;

while norm(x, inf) > 1
    % adjust the near spt terms
    for i = 2:bitplaneDepth
        if x(i-1) > 0.0 && x(i) < 0.0
            x(i-1) = x(i-1) - 1;
            x(i) = x(i) + 2;
        elseif x(i-1) < 0.0 && x(i) > 0.0
            x(i-1) = x(i-1) + 1;
            x(i) = x(i) - 2;
        end
    end
    % adjust for absolute values greater than 1
    for i = bitplaneDepth:-1:2
        x(i-1) = x(i-1) + sign(x(i))*floor(abs(x(i))/2);
        x(i) = x(i) - 2*sign(x(i))*floor(abs(x(i))/2);
    end
    % readjust the near spt terms
    for i = 2:bitplaneDepth
        if x(i-1) > 0.0 && x(i) < 0.0
            x(i-1) = x(i-1) - 1;
            x(i) = x(i) + 2;
        elseif x(i-1) < 0.0 && x(i) > 0.0
            x(i-1) = x(i-1) + 1;
            x(i) = x(i) - 2;
        end
    end
    % overflow
    if (abs(x(1)) >= 2 && count == 1) || (abs(x(1)) > 2)
        x = sign(x(1))*ones(1, bitplaneDepth);
        fprintf('REAPPROX_CONST: overflow.\n')
    elseif abs(x(1)) == 2 && count == 0
        x_aux = sign(x(1))*ones(1, bitplaneDepth);
        x(1) = x(1) - sign(x(1));
        x = x + x_aux;
        count = 1;
    end
end
nz = nnz(x);
while nz > maxSumsfilt  % if the number of non-zeros is greater than maxSums...
    for it = bitplaneDepth:-1:1
        if x(it) ~= 0
            x(it) = 0;
            nz = nnz(x);
            break
        end
    end
end
y = x;



% EoF

