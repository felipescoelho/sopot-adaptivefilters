function [y_prio, e_prio, coefficients,...
    y_depth, e_depth, desired_depth, w_depth] = RLS_spt(desired, input, S)

%   RLS_spt.m
%       Implements the Alternative RLS Algorithm for an approximation in
%       Sums of Powers of Two -- SPT. The RLS_Alt algorithm uses an
%       auxiliar variable (psi) in order to reduce the computational
%       burden. This proposed SPT variant uses a SPT indexing in order to
%       reduce the computational complexity.
%       (The reference RLS_Alt algorithm comes from, Algorithm 5.4 - book:
%       Adaptive Filtering: Algorithms and Practical Implementation, Diniz)
%
%
%   Syntax:
%       [y_priori, e_priori, coefficients] = RLS_spt(desired, input, S)
%
%
%   Input Arguments:
%       . desired   : Desired signal.
%                                                   COLUMN(nIterations, 1)
%       . input     : Input signal.
%                                                   COLUMN(nIterations, 1)
%       . S         : Structure with the following fields:
%           - filterOrder   : Order of the FIR filter.
%           - initCoeff  : Initial filter coefficients.
%                                               COLUMN(nCoefficients, 1)
%           - delta         : The matrix delta*eye is the initial value of
%                             the inverse of the deterministic
%                             autocorrelation matrix.
%           - lambda        : Forgetting factor. (0 << lambda < 1)
%           - bitplaneDepth : Maximun bit depth for the spt approximation.
%           - maxSums       : Maximum number of non-zero terms in the spt
%                             approximation for the input signal.
%           - maxSumsfilt   : Maximum number of non-zero terms in the spt
%                             approximation for the adaptive filter
%                             operations.
%           - maxPot        : Greater power for the spt approximation
%
%
%   Output Arguments:
%       . y_priori      : Estimated output for each iteration.
%                                                   COLUMN(nIterations, 1)
%       . e_priori      : Error for each iteration.
%                                                   COLUMN(nIterations, 1)
%       . coefficients  : Estimated coefficients for each iteration.
%                                   2-D ARRAY(nCoefficients, nIterations)
%
%
%   Author: Luiz Felipe da S. Coelho - luizfscoelho.92@gmail.com
%
%


%   Initialization:
bitplaneDepth = S.wordLength;
maxSums = S.maxSums;
maxSumsfilt = S.maxSumsFilt;
maxPot = S.maxPot;
nCoefficients = S.filterOrderNo + 1;
nIterations = length(input);


%   SPT approximation:
desired_spt = spt_approx(desired, bitplaneDepth, maxSums, maxPot);
initCoeff_spt = spt_approx(S.initCoeffs, bitplaneDepth, maxSumsfilt, maxPot);
delta_spt = spt_approx(S.delta, bitplaneDepth, maxSumsfilt, maxPot);
lambda_spt = spt_approx(S.lambda, bitplaneDepth, maxSumsfilt, maxPot);


%   Pre-allocation:
desired_depth = zeros(nIterations, 1);
%       a priori
y_prio_spt = zeros(nIterations, bitplaneDepth);
y_depth = zeros(nIterations, 1);
e_prio_spt = zeros(nIterations, bitplaneDepth);
e_depth = zeros(nIterations, 1);
%       a posteriori
y_post_spt = zeros(nIterations, bitplaneDepth);
e_post_spt = zeros(nIterations, bitplaneDepth);
%       coefficients
coefficients_spt = zeros(nCoefficients, bitplaneDepth, (nIterations+1));
w_depth = -200*ones(nIterations, nCoefficients);
%       S_d
S_d = zeros(nCoefficients, nCoefficients, bitplaneDepth);


%   Initial state:
%       coefficients
coefficients_spt(:, :, 1) = initCoeff_spt;
%       S_d(i-1)
for j = 1:nCoefficients
    S_d(j, j, :) = delta_spt;
end


%   Improve sorce code regularity:
extendedInput = [zeros(nCoefficients-1, 1)
                 input];


%   Body:
for i = 1:nIterations
    %   tapped delay line
    regressor = extendedInput(i+nCoefficients-1:-1:i, :);
    regressor_spt = spt_approx(regressor, bitplaneDepth, maxSums, maxPot);
    %   a priori output estiomation
    y_prio_spt(i, :) = p_rowcol(regressor_spt.',...
        coefficients_spt(:, :, i), maxSumsfilt, maxPot);
    for m = 0:bitplaneDepth-1
        y_depth(i) = -200;
        if y_prio_spt(i, bitplaneDepth - m) ~= 0
            y_depth(i) = maxPot - (bitplaneDepth-m) + 1;
            break
        end
    end
    %   a prior error estimation
    e_prio_aux = desired_spt(i, :) - y_prio_spt(i, :);
    for m = 0:bitplaneDepth-1
        desired_depth(i) = -200;
        if desired_spt(i, bitplaneDepth - m) ~= 0
            desired_depth(i) = maxPot - (bitplaneDepth-m) + 1;
            break
        end
    end
    e_prio_spt(i, :) = reapprox_const(e_prio_aux, maxSumsfilt);
    for m = 0:bitplaneDepth-1
        w_depth(i) = -200;
        if e_prio_spt(i, bitplaneDepth - m) ~= 0
            e_depth(i) = maxPot - (bitplaneDepth-m) + 1;
            break
        end
    end
    %   psi
    psi = p_matcol(S_d, regressor_spt, maxSumsfilt, maxPot);
    %   S_d(k)
    num = p_colrow(psi, psi.', maxSumsfilt, maxPot);
    den_aux1 = lambda_spt + p_rowcol(psi.', regressor_spt, maxSumsfilt,...
        maxPot);
    den_aux2 = reapprox_const(den_aux1, maxSumsfilt);
    den = invert_spt(den_aux2, maxSumsfilt, maxPot);  % must be in the format [1, bitplaneDepth]
    S_d_aux1 = p_matconst(num, den, maxSumsfilt, maxPot);
    S_d_aux2 = S_d - S_d_aux1;
    S_d_aux3 = reapprox_mat(S_d_aux2, maxSumsfilt);
    adbmal = invert_spt(lambda_spt, maxSumsfilt, maxPot);
    S_d = p_matconst(S_d_aux3, adbmal, maxSumsfilt, maxPot);
    %   coefficient adjustment
    coefficients_aux1 = p_matconst(S_d, e_prio_spt(i, :), maxSumsfilt,...
        maxPot);
    coefficients_aux2 = p_matcol(coefficients_aux1, regressor_spt,...
        maxSumsfilt, maxPot);
    coefficients_aux3 = coefficients_spt(:, :, i) + coefficients_aux2;
    coefficients_spt(:, :, i+1) = reapprox_col(coefficients_aux3,...
        maxSumsfilt);
    for coefs = 1:nCoefficients
        for m = 0:bitplaneDepth-1
            if coefficients_spt(coefs, bitplaneDepth - m, i+1) ~= 0
                w_depth(i, coefs) = maxPot - (bitplaneDepth-m) + 1;
                break
            end
        end
    end
    %   a posteriori output estimation
    y_post_spt(i, :) = p_rowcol(coefficients_spt(:, :, i+1).',...
        regressor_spt, maxSumsfilt, maxPot);
    e_post_aux = desired_spt(i, :) - y_post_spt(i, :);
    e_post_spt(i, :) = reapprox_const(e_post_aux, maxSumsfilt);
end


%   Recovering output:
%       signals
y_prio = spt_signal_recover(y_prio_spt, maxPot);
e_prio = spt_signal_recover(e_prio_spt, maxPot);
%       coefficients
coefficients = spt_coefficient_recover(coefficients_spt, maxPot);


% EoF


