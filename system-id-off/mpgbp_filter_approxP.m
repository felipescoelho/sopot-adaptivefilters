function [hout, km] = mpgbp_filter_approxP(hin, nAdders, P)

% 
% This function acompains the mansucript "FIR Filter Design Based on Successive
% Approximation of Vectors" by Da Silva, Lovisolo, Dutra and Diniz.
%
% function hout = mpgbp_filter_approxP(hin, nAdders, P) ;
%
% It approximates the filter hin for a limmiting bit-depth or POT terms given
% by Bmax, using a dictionary defined with P non-zero elements (+1,-1).
%
% *** Output:
%
% hout  -> filter approximation using the above limitants
%
% *** Input:
%
% hin      - original filter design
% nAdders  - max number of POT terms to be used in the approximation
% P        - number of non-zero codeword components
%
% *** Usage example:
% 
% Being hin a FIR filter use
%
% hout = mpgbp_filter_approxP(hin, nAdders, P)
%
% To obtain the approximation of filter hin using SOPOTs for a maximum depth 
% of POT terms nAdders, using the dictionary defined by P.
%
if nnz(hin) == 0
    hout = hin';
    km = 0;
else

[~, ncolumn] = size(hin) ;

if ncolumn == 1
	hin = hin';
end

[~, ncolumn] = size(hin) ;
N = ncolumn ;				  % Number of filter taps (= vector dimension)

hout = zeros(1,N) ;

km = 0;

% Codebook setup:
if P >= N
    return
end

% P = floor(sqrt(N)) ;	% Number of non-zero components in each codeword

h_residue   = hin ;
norm(h_residue) ;

% MPGBP approximation process:
m_idx       = zeros(1,P) ;          % Array containing indices of P largest components of hin (in absolute value)
alpha_R     = 0.5 ;                 % Initial reproduction threshold
error       = norm(hin) ;           % Initial reproduction error

count       = 0 ;                        % current number of POT terms used in the approx.
comp_count  = zeros(1,N) ;          % Number of POT terms per component

hpast       = zeros(1,N) ;              % Output filter buffer (keeps a copy of previous pass approximation)

contador    = 0;
while count < nAdders
    v_rm      = zeros(1, N) ;                        % Approximation codeword
    temp_h    = h_residue ;

    %   Find reproduction vector:
    [~, idx_temp]    = sort(abs(temp_h),'descend') ;
    m_idx            = idx_temp(1:P) ;
    for i=1:P
        v_rm (m_idx(i)) = 1.0 * sign(temp_h(m_idx(i))) ;
    end

	v_rn         = v_rm./norm(v_rm) ;
    ip           = (h_residue * v_rn') / (norm(v_rm)) ;
    

    km           = -ceil(log2(3.0/(4.0*ip))) ;
    kfactor      = 2.0 ^ km ;

    h_residue    = h_residue - kfactor * v_rm ;
    hout         = hout + kfactor * v_rm ;

    [count, ~]   = adder_count(hout) ;
%     count
      contador = contador +1;
    if count >= nAdders
        hout = hpast ;
        break
    else
        hpast = hout ;
    end

    if norm(h_residue) < 1e-10
        break ;
    end

%     pause
end
% h_residue = hin - hout ;
% error = norm(h_residue) ;

end
end
