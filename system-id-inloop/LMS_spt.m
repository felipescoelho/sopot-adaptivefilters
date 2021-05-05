function [y, e, w] = LMS_spt(d, x, S)

% LMS_spt
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
mu = S.step;
K = length(d);


% SOPOT approximation
d_spt = spt_approx(d, word_len, max_sum, max_pot);
w_i_spt = spt_approx(S.initCoeffs, word_len, max_sum, max_pot);
mu_spt = spt_approx(mu, word_len, max_sum, max_pot);


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
    step_aux = inner_p(mu_spt, e_spt(k, :), maxs_filt, max_pot);
    w_aux = w_spt(:, :, k) +  p_constcol(step_aux, regressor_spt, maxs_filt, max_pot);
    w_spt(:, :, k+1) = reapprox_col(w_aux, maxs_filt);
end


% recover
e = spt_signal_recover(e_spt, max_pot);
y = spt_signal_recover(y_spt, max_pot);
w = spt_coefficient_recover(w_spt, max_pot);


% EoF


