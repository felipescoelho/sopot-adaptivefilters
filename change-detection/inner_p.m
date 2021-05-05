function z = inner_p(x, y, maxSumsfilt, maxPot)

%   inner_p.m
%       Implements the multiplication of two numbers approximated by the
%       spt, this operation is a key-operation for the implementation of
%       the spt-RLS.
%
%
%   Syntax:
%       z = inner_p(x, y)
%
%
%   Input Arguments:
%       . x     : array containing the spt indexes.
%                                                   ROW(1, bitplaneDepth)
%       . y     : array containing the spt indexes.
%                                                   ROW(1, bitplaneDepth)
%
%
%   Output Arguments:
%       . z     : array contaning the spt indexes.
%                                                   ROW(1, bitplaneDepth)
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%

[~, bitplaneDepth] = size(x);

%   Computing Operation:

z_aux1 = zeros(bitplaneDepth, 2*bitplaneDepth - 1);
for i = 1:bitplaneDepth
    for j = 1:bitplaneDepth
        z_aux1(i, i+j-1) = x(i)*y(j);
    end
end


%   Reapproximation Process:

z_aux2 = sum(z_aux1, 1);
count = 0;
clear z_aux1


while norm(z_aux2, inf) > 1.0
    % adjusting near spt terms
    for i = 2:(2*bitplaneDepth - 1)
        if z_aux2(i-1) > 0.0 && z_aux2(i) < 0.0
            z_aux2(i-1) = z_aux2(i-1) - 1;
            z_aux2(i) = z_aux2(i) + 2;
        elseif z_aux2(i-1) < 0.0 && z_aux2(i) > 0.0
            z_aux2(i-1) = z_aux2(i-1) + 1;
            z_aux2(i) = z_aux2(i) - 2;
        end
    end
    % adjusting for absolute values grater than 1
    for i = (2*bitplaneDepth - 1):-1:2
        z_aux2(i-1) = z_aux2(i-1) + sign(z_aux2(i))*floor(abs(z_aux2(i))/2);
        z_aux2(i) = z_aux2(i) - 2*sign(z_aux2(i))*floor(abs(z_aux2(i))/2);
    end
    % adjusting near spt terms
    for i = 2:(2*bitplaneDepth-1)
        if z_aux2(i-1) > 0.0 && z_aux2(i) < 0.0
            z_aux2(i-1) = z_aux2(i-1) - 1;
            z_aux2(i) = z_aux2(i) + 2;
        elseif z_aux2(i-1) < 0.0 && z_aux2(i) > 0.0
            z_aux2(i-1) = z_aux2(i-1) + 1;
            z_aux2(i) = z_aux2(i) - 2;
        end
    end
    % overflow
    my_ovf = z_aux2(1:maxPot);
    if nnz(my_ovf) ~= 0
        if abs(my_ovf(end)) == 1 && count == 0
            z_aux3 = sign(z_aux2(maxPot))*ones(1, 2*bitplaneDepth-1-maxPot);
            z_aux2(maxPot) = z_aux2(maxPot) - sign(z_aux2(maxPot));
            z_aux4 = zeros(1, 2*bitplaneDepth - 1);
            z_aux4(maxPot+1:(2*bitplaneDepth - 1)) = z_aux3;
            z_aux2 = z_aux2 + z_aux4;
            count = 1;
        else
            %fprintf('INNER_P: overflow.\n')
            for i = 1:maxPot
                if my_ovf(i) ~= 0
                    z_aux2 = sign(my_ovf(i))*ones(1, (2*bitplaneDepth - 1));
                    break
                end
            end
        end
    end
end

z = z_aux2((maxPot+1):(bitplaneDepth+maxPot));
nz = nnz(z);

%  Imposing the complexity limit.
while nz > maxSumsfilt
    for it = bitplaneDepth:-1:1
        if z(it) ~= 0
            z(it) = 0;
            nz = nnz(z);
            break
        end
    end
end
% roundoff
if nz == 0
    %fprintf('INNER_P: roundoff.\n')
    for i = (bitplaneDepth+maxPot+1):(2*bitplaneDepth - 1)
        if z_aux2(i) ~= 0
            z(end) = z_aux2(i);
            break
        end
    end
end


% EoF



