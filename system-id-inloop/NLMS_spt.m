function [y, e, w] = NLMS_spt(d, x, S)

% NLMS_spt
%   Implements the LMS algorithm using SOPOT representation.
%
%
%   Author: Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
%


% initialization
word_len = S.wordLength;
max_sum = S.maxSums;
maxs_filt = S.maxSumsFilt;
max_pot = S.maxPot;
nCoeff = S.filterOrderNo + 1;
K = length(d);


% SOPOT approximation
d_spt = spt_approx(d, word_len, max_sum, max_pot);
w_i_spt = spt_approx(S.initCoeffs, word_len, max_sum, max_pot);
mu = spt_approx(S.step, word_len, max_sum, max_pot);
gamma = spt_approx(S.gamma, word_len, max_sum, max_pot);


% memory allocation
e_spt = zeros(K, word_len);
y_spt = zeros(K, word_len);
w_spt = zeros(nCoeff, word_len, K+1);


% initial state
w_spt(:, :, 1) = w_i_spt;


% prefix input
prefixedInput = [zeros(nCoeff-1, 1)
                 x];


% computing
for k = 1:K
    % tapped delay line
    regressor = prefixedInput(k+nCoeff-1:-1:k, :);
    regressor_spt = spt_approx(regressor, word_len, max_sum, max_pot);
    % output signal
    y_spt(k, :) = p_rowcol(regressor_spt', w_spt(:, :, k), maxs_filt, max_pot);
    % error signal
    e_aux = d_spt(k, :) - y_spt(k, :);
    e_spt(k, :) = reapprox_const(e_aux, maxs_filt);
    % coefficient update
    num_aux = inner_p(mu, e_spt(k, :), maxs_filt, max_pot);
    num = p_constcol(num_aux, regressor_spt, maxs_filt, max_pot);
    den_aux = p_rowcol(regressor_spt', regressor_spt, maxs_filt, max_pot);
    den = reapprox_const(gamma + den_aux, maxs_filt);
    ned = invert_spt(den, maxs_filt, max_pot);
    w_aux = w_spt(:, :, k) + p_constcol(ned, num, maxs_filt, max_pot);
    w_spt(:, :, k+1) = reapprox_col(w_aux, maxs_filt);
end


% recover
e = spt_signal_recover(e_spt, max_pot);
y = spt_signal_recover(y_spt, max_pot);
w = spt_coefficient_recover(w_spt, max_pot);


% EoF


