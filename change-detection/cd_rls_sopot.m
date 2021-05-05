function [alarm_bool] = cd_rls_sopot(d, maxSums, maxPot, maxSumsFilt,...
                                     threshold, S)

%   cd_rls_sopot.m
%       Implements a change detector based on a RLS adaptive filter using
%       sopot representation.
%
%   Syntax:
%       [alarm] = cd_rls_sopot(y, S);
%
%   Input args.:
%       . y     : signal to be evaluated. (COLUMN)
%       . S     : structure with the following fields
%           - filterOrderNo         : order of the FIR filter.
%           - initialCoefficients   : initial filter coefficients (COLUMN)
%           - lambda                : forgetting factor. (0 << lambda < 1)
%           - delta                 : initial value of the inverse of the
%                                     deterministic autocorrelation matrix.
%           - windowSize            : window size for the mean filter.
%           - threshold             : threshold for the alarm.
%           - bitplaneDepth         : maximum bit depth for the SOPOT
%                                     approximation.
%           - maxSums               : maximum number of non-zero terms in
%                                     the SOPOT approximation for the input
%                                     signal.
%           - maxSumsFilt           : maximum number of non-zero terms in
%                                     the SOPOT approximation for the
%                                     adaptive filter operations.
%           - maxPot                : greatest power for the SOPOT
%                                     approximation.
%
%   Outpu args.:
%       . alarm : binary mask with detection through time.
%
%
%   Author:
%       . Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
%
%


% initialization
word_len = S.bitplaneDepth;
active_bits = maxSums;
internal_active_bits = maxSumsFilt;
max_pot = maxPot;
num_coeff = S.filterOrderNo + 1;
win_size = S.windowSize;
num_iter = length(d);


% sopot approximation
d_sopot = spt_approx(d, word_len, active_bits, max_pot);
h_sopot = spt_approx(threshold, word_len, active_bits, max_pot);
h = spt_signal_recover(h_sopot, max_pot);
delta_sopot = spt_approx(S.delta, word_len, active_bits, max_pot);
lambda_sopot = spt_approx(S.lambda, word_len, active_bits, max_pot);
ezis_niw = spt_approx(1/win_size, word_len, internal_active_bits, max_pot);
i_coeff_sopot = spt_approx(S.initialCoefficients, word_len, active_bits,...
                           max_pot);

% pre allocation
epsilon_sopot = zeros(num_iter, word_len);
d_hat = zeros(num_iter, word_len);
alarm = ones(1, num_iter);
coeff_vector = zeros(num_coeff, word_len, num_iter+1);
S_d = zeros(num_coeff, num_coeff, word_len);
adbmal = invert_spt(lambda_sopot, internal_active_bits, max_pot);



% initial state
coeff_vector(:, :, 1) = i_coeff_sopot;
for j = 1:num_coeff
    S_d(j, j, :) = delta_sopot;
end


% delay input
x = zeros(num_iter, 1);
x(2:end, :) = d(1:end-1, :);


% adjusting regularity
extendedInput = [zeros(num_coeff-1, 1)
                 x];
extendedEpsilon = [zeros(win_size-1, word_len)
                   epsilon_sopot];

% body
for i = 1:num_iter
    % delay line
    regressor = extendedInput(i+num_coeff-1:-1:i);
    regressor_sopot = spt_approx(regressor, word_len, active_bits, max_pot);
    
    % output -> y = d_hat
    d_hat(i, :) = p_rowcol(regressor_sopot.', coeff_vector(:, :, i),...
                           internal_active_bits, max_pot);
    
    % prediction error
    epsilon_aux = d_sopot(i, :) - d_hat(i, :);
    extendedEpsilon(i+win_size-1, :) = reapprox_const(epsilon_aux,...
                                                    internal_active_bits);
    
    % prediction error evaluation
    epsilon_regressor = extendedEpsilon(i+win_size-1:-1:i, :);
    epsilon_regressor_aux = p_rowconst(epsilon_regressor.', ezis_niw,...
                                       internal_active_bits, max_pot);
    sqrd_epsilon_sopot = p_rowcol(epsilon_regressor.',...
                                  epsilon_regressor_aux.',...
                                  internal_active_bits, max_pot);
    sqrd_epsilon = spt_signal_recover(sqrd_epsilon_sopot, max_pot);
    if sqrd_epsilon >= h
        alarm(i) = 0;
    end
    
    % psi
    psi = p_matcol(S_d, regressor_sopot, internal_active_bits, max_pot);
    
    % inverse of the deterministic autocorrelation matrix - S_d
    num = p_colrow(psi, psi.', internal_active_bits, max_pot);
    den_aux1 = lambda_sopot + p_rowcol(psi.', regressor_sopot,...
                                       internal_active_bits, max_pot);
    den_aux2 = reapprox_const(den_aux1, internal_active_bits);
    den = invert_spt(den_aux2, internal_active_bits, max_pot);
    S_d_aux1 = p_matconst(num, den, internal_active_bits, max_pot);
    S_d_aux2 = S_d - S_d_aux1;
    S_d_aux3 = reapprox_mat(S_d_aux2, internal_active_bits);
    S_d = p_matconst(S_d_aux3, adbmal, internal_active_bits, max_pot);
    
    % coefficients
    coeff_aux1 = p_matconst(S_d, extendedEpsilon(i+win_size-1, :),...
                            internal_active_bits, max_pot);
    coeff_aux2 = p_matcol(coeff_aux1, regressor_sopot,...
                          internal_active_bits, max_pot);
    coeff_aux3 = coeff_vector(:, :, i) + coeff_aux2;
    coeff_vector(: ,:, i+1) = reapprox_col(coeff_aux3, internal_active_bits);
end

alarm_bool = alarm == 1;

% EoF


