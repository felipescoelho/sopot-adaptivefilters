function y = spt_approx(x, bitplaneDepth, maxSums, maxPot)

%   spt_approx.m
%       This function approximates the signal x in a limiting bit plane
%       depth and a limiting power of two terms, using a dictionary with P
%       non-zero elements (-1, 1). Based on Lisandro's
%       MPGBP_filter_approxP.m
%
%
%   Syntax:
%       y = spt_approx(x, bitplaneDepth, maxSums)
%
%
%   Input Arguments:
%       . x                 : Original signal
%       . bitplaneDepth     : Maximum depth for the bit-plane.
%       . maxSums           : Maximum number of sums for the spt
%                             approximation.
%
%
%   Output Arguments:
%       . y     : Approximated signal in the spt format.
%                                   2-D ARRAY(signalLength, bitplaneDepth)
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%


len = length(x);
P = floor(sqrt(len));
y = zeros(len, bitplaneDepth);
x_residue = x;
x_spt = zeros(len, 1);

if nnz(x) == 0
    return
end

while norm(x_residue) > 1e-12
    v_rm = zeros(len, 1);  % approximation codeword
    temp_x = x_residue;
    
    
    %   Finding reproduction vector
    [~, idx_temp] = sort(abs(temp_x), 'descend');
    m_idx = idx_temp(1:P);
    
    for i = 1:P
        v_rm(m_idx(i)) = 1.0*sign(temp_x(m_idx(i)));
    end
    
    norm_v_rm = norm(v_rm);
    v_rn = v_rm./norm_v_rm;
    ip = (x_residue'*v_rn)./norm_v_rm;
    km = ceil(log2(3.0/(4.0*ip)));
    
    kfactor = 2.0^-km;
    x_residue = x_residue - kfactor*v_rm;
    x_spt = x_spt + kfactor*v_rm;
    
    %   spt index allocation:
    if (km+maxPot+1) <= bitplaneDepth && km >= -maxPot
        for i = 1:P
            y(m_idx(i), km+maxPot+1) = y(m_idx(i), km+maxPot+1) + v_rm(m_idx(i));
        end
    end
end

y = reapprox_col(y, maxSums);

return



% EoF

