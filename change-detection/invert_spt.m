function y = invert_spt(x, maxSumsfilt, maxPot)

%   invert_spt.m
%       This function does an approximation of the inversion of a spt
%       approximated number. It uses the Maclaurin Series for the
%       approximation in the case which the number is within its
%       constraints.
%
%
%   Syntax:
%       y = invert_spt(x)
%
%
%   Input Arguments:
%       . x     : number to be inverted.
%                                                   ROW (1, bitplaneDepth)
%
%
%   Output Arguments:
%       . y     : inverted number.
%                                                   ROW (1, bitplaneDepth)
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%

bitplaneDepth = length(x);
if maxPot == 1.0
    if x(1) == 1.0
        y = zeros(1, bitplaneDepth);
        y(3) = 1;
       % fprintf('The inverted number is greater than 2\n')
    elseif x(1) == -1
        y = zeros(1, bitplaneDepth);
        y(3) = -1;
      %  fprintf('The inverted number is smaller than -2\n')
    else
        um = zeros(1, bitplaneDepth);
        um(2) = 1;
        x_2 = inner_p((um - x),(um - x), maxSumsfilt, maxPot);
        x_3 = inner_p(x_2, (um - x), maxSumsfilt, maxPot);
        y = um + um - x + x_2 + x_3;
        y = reapprox_const(y, maxSumsfilt);
    end
elseif maxPot == 2.0
    if x(1) == 1.0
        y = zeros(1, bitplaneDepth);
        y(5) = 1;
       % fprintf('The inverted number is greater than 4\n')
    elseif x(1) == -1
        y = zeros(1, bitplaneDepth);
        y(5) = -1;
      %  fprintf('The inverted number is smaller than -4\n')
    elseif x(2) == 1.0
        y = zeros(1, bitplaneDepth);
        y(4) = 1;
        %fprintf('The inverted number is greater than 2\n')
    elseif x(2) == -1
        y = zeros(1, bitplaneDepth);
        y(4) = -1;
        %fprintf('The inverted number is smaller than -2\n')
    else
        um = zeros(1, bitplaneDepth);
        um(3) = 1;
        x_2 = inner_p((um - x),(um - x), maxSumsfilt, maxPot);
        x_3 = inner_p(x_2, (um - x), maxSumsfilt, maxPot);
        y = um + um - x + x_2 + x_3;
        y = reapprox_const(y, maxSumsfilt);
    end
elseif maxPot == 3.0
    if x(1) == 1.0
        y = zeros(1, bitplaneDepth);
        y(7) = 1;
        %fprintf('The inverted number is greater than 8\n')
    elseif x(1) == -1
        y = zeros(1, bitplaneDepth);
        y(7) = -1;
        %fprintf('The inverted number is smaller than -8\n')
    elseif x(2) == 1
        y = zeros(1, bitplaneDepth);
        y(6) = 1;
        %fprintf('The inverted number is greater than 4\n')
    elseif x(2) == -1
        y = zeros(1, bitplaneDepth);
        y(6) = -1;
        %fprintf('The inverted number is smaller than -4\n')
    elseif x(3) == 1
        y = zeros(1, bitplaneDepth);
        y(5) = 1;
        %fprintf('The inverted number is greater than 2\n')
    elseif x(3) == -1
        y = zeros(1, bitplaneDepth);
        y(5) = -1;
        %fprintf('The inverted number is smaller than -2\n')
    else
        um = zeros(1, bitplaneDepth);
        um(4) = 1;
        x_2 = inner_p((um - x),(um - x), maxSumsfilt, maxPot);
        x_3 = inner_p(x_2, (um - x), maxSumsfilt, maxPot);
        y = um + um - x + x_2 + x_3;
        y = reapprox_const(y, maxSumsfilt);
    end
elseif maxPot == 4.0
    if x(1) == 1.0
        y = zeros(1, bitplaneDepth);
        if bitplaneDepth > 8
            y(9) = 1;
        end
        %fprintf('The inverted number is greater than 16\n')
    elseif x(1) == -1.0
        y = zeros(1, bitplaneDepth);
        if bitplaneDepth > 8
            y(9) = -1;
        end
        %fprintf('The inverted number is smaller than -16\n')
    elseif x(2) == 1.0
        y = zeros(1, bitplaneDepth);
        y(8) = 1;
        %fprintf('The inverted number is greater than 8\n')
    elseif x(2) == -1
        y = zeros(1, bitplaneDepth);
        y(8) = -1;
        %fprintf('The inverted number is smaller than -8\n')
    elseif x(3) == 1
        y = zeros(1, bitplaneDepth);
        y(7) = 1;
        %fprintf('The inverted number is greater than 4\n')
    elseif x(3) == -1
        y = zeros(1, bitplaneDepth);
        y(7) = -1;
        %fprintf('The inverted number is smaller than -4\n')
    elseif x(4) == 1
        y = zeros(1, bitplaneDepth);
        y(6) = 1;
        %fprintf('The inverted number is greater than 2\n')
    elseif x(4) == -1
        y = zeros(1, bitplaneDepth);
        y(6) = -1;
        %fprintf('The inverted number is smaller than -2\n')
    else
        um = zeros(1, bitplaneDepth);
        um(5) = 1;
        x_2 = inner_p((um - x),(um - x), maxSumsfilt, maxPot);
        x_3 = inner_p(x_2, (um - x), maxSumsfilt, maxPot);
        y = um + um - x + x_2 + x_3;
        y = reapprox_const(y, maxSumsfilt);
    end
elseif maxPot == 5
    if x(1) == 1.0
        y = zeros(1, bitplaneDepth);
        y(11) = 1;
        %fprintf('The inverted number is greater than 32\n')
    elseif x(1) == -1.0
        y = zeros(1, bitplaneDepth);
        y(11) = -1;
        %fprintf('The inverted number is smaller than -32\n')
    elseif x(2) == 1.0
        y = zeros(1, bitplaneDepth);
        y(10) = 1;
        %fprintf('The inverted number is greater than 16\n')
    elseif x(2) == -1
        y = zeros(1, bitplaneDepth);
        y(10) = -1;
        %fprintf('The inverted number is smaller than -16\n')
    elseif x(3) == 1
        y = zeros(1, bitplaneDepth);
        y(9) = 1;
        %fprintf('The inverted number is greater than 8\n')
    elseif x(3) == -1
        y = zeros(1, bitplaneDepth);
        y(9) = -1;
        %fprintf('The inverted number is smaller than -8\n')
    elseif x(4) == 1
        y = zeros(1, bitplaneDepth);
        y(8) = 1;
        %fprintf('The inverted number is greater than 4\n')
    elseif x(4) == -1
        y = zeros(1, bitplaneDepth);
        y(8) = -1;
        %fprintf('The inverted number is smaller than -4\n')
    elseif x(5) == 1
        y = zeros(1, bitplaneDepth);
        y(7) = 1;
        %fprintf('The inverted number is greater than 2\n')
    elseif x(5) == -1.0
        y = zeros(1, bitplaneDepth);
        y(7) = -1;
        %fprintf('The inverted number is smaller than -2\n')
    else
        um = zeros(1, bitplaneDepth);
        um(6) = 1;
        x_2 = inner_p((um - x),(um - x), maxSumsfilt, maxPot);
        x_3 = inner_p(x_2, (um - x), maxSumsfilt, maxPot);
        y = um + um - x + x_2 + x_3;
        y = reapprox_const(y, maxSumsfilt);
    end
end



% EoF

