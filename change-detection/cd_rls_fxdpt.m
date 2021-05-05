function [alarm_bool] = cd_rls_fxdpt(y, wordLength, fractionLength,...
                                     threshold, S)

%   cd_rls_fxdpt.m
%       Implements a change detector based on a RLS adaptive filter using
%       fixed-point representation.
%
%   Syntax:
%       [alarm] = cd_rls_fxdpt(y, S);
%
%   Input args.:
%       . y     : signal to be evaluated. (COLUMN)
%       . S     : structure with the following fields
%           - filterOrderNo         : order of the FIR filter.
%           - initialCoefficients   : initial filter coefficients (COLUMN)
%           - delta                 : the matrix delta*eye is the initial
%                                     value of the inverse of the
%                                     deterministic autocorrelation matrix.
%           - lambda                : forgetting factor. (0 << lambda < 1)
%           - windowSize            : window size for the mean filter.
%           - threshold             : threshold for the alarm.
%           - wordLength            : wordlength for the fixed-point.
%           - fractionLength        : number of bits after the binary point
%
%   Output args.:
%       . alarm : binary mask with detection through time.
%
%
%   Author:
%       . Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
%
%


% initialization
num_coeff = S.filterOrderNo+1;
num_iter = length(y);
word_len = wordLength;
frac_len = fractionLength;
win_size = S.windowSize;


% fixed-point
divisor = fi(1/win_size, 1, word_len, frac_len);
d = fi(y, 1, word_len, frac_len);
h = fi(threshold, 1, word_len, frac_len);
lambda = fi(S.lambda, 1, word_len, frac_len);

% pre allocation
epsilon = fi(zeros(num_iter, 1), 1, word_len, frac_len);
coeff_vector = fi(zeros(num_coeff, num_iter+1), 1, word_len, frac_len);
d_hat = fi(zeros(num_iter, 1), 1, word_len, frac_len);
alarm = ones(1, num_iter);


% initial state
coeff_vector(:, 1) = fi(S.initialCoefficients, 1, word_len, frac_len);
S_d = fi(S.delta*eye(num_coeff), 1, word_len, frac_len);


% delay input
x = fi(zeros(num_iter, 1), 1, word_len, frac_len);
x(2:end) = d(1:end-1);


% regularity
extendedInput = fi([zeros(num_coeff-1, 1)
                    x], 1, word_len, frac_len);
extendedEpsilon = fi([zeros(win_size-1, 1)
                      epsilon], 1, word_len, frac_len);


% body
for i = 1:num_iter
    % delay line
    regressor = extendedInput(i+(num_coeff-1):-1:i);
    % estimation for d -> d_hat
    d_hat(i, :) = coeff_vector(:, i)'*regressor;
    % prediction error
    extendedEpsilon(i+win_size-1) = d(i) - d_hat(i);
    % prediction error eval
    epsilon_regressor = extendedEpsilon(i+win_size-1:-1:i);
    sqrd_epsilon = epsilon_regressor.'*epsilon_regressor;
    mean_epsilon = sqrd_epsilon*divisor;
    if mean_epsilon >= h
        alarm(i) = 0;
    end
    % psi
    psi = S_d*regressor;
    % inverse of the deterministic autocorrelation matrix - S_d
    S_d(:) = (1/lambda)*(S_d - (psi*psi')/(lambda + psi'*regressor));
    % coefficients
    coeff_vector(:, i+1) = coeff_vector(:, i) +...
        extendedEpsilon(i+win_size-1)*S_d*regressor;
end

alarm_bool = alarm == 1;


% EoF


