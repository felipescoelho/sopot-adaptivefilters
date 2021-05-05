function y = spt_signal_recover(x, maxPot)

%   spt_signal_recover.m
%       This function recovers a signal in the spt approximation format
%       back into the conventional numeric representation.
%
%
%   Syntax:
%       y = spt_signal_recover(x)
%
%
%   Output:
%       . y     : Signal recovered from the spt approximation.
%                                                   COLUMN(signalLength, 1)
%
%
%   Input:
%       . x     : Signal to be recovered.
%                                   2-D ARRAY(signalLength, bitplaneDepth)
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%


[len, bitplaneDepth] = size(x);
y = zeros(len, 1);

for it = 1:len
    y_aux = zeros(1, bitplaneDepth);
    for k = 1:bitplaneDepth
        y_aux(1, k) = x(it, k)*2^((maxPot+1) - k);
        %{
        if x(it, k) == 1.0
            y_aux(1, k) = 2^(-k+maxPot+1);
        elseif x(it, k) == -1.0
            y_aux(1, k) = -2^(-k+maxPot+1);
        end
        %}
    end
    y(it) = sum(y_aux, 2);
end

return




% EoF

