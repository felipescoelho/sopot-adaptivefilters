function    [y, e, w] =   NLMS_fixed_point(d, x, S)

%   NLMS.m
%       Implements the Normalized LMS algorithm for COMPLEX valued data.
%       (Algorithm 4.3 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
%
% - Update: 25.05.2020 (dd.mm.yyyy)
%   Modified by Luiz Felipe da S. Coelho to be used as a fixed-point LMS
%                                    - email: luizfelipe.coelho@smt.ufrj.br
%
%   Syntax:
%       [outputVector,errorVector,coefficientVector] = NLMS(desired,input,S)
%
%   Input Arguments:
%       . desired   : Desired signal.                               (ROW vector)
%       . input     : Signal fed into the adaptive filter.          (ROW vector)
%       . S         : Structure with the following fields
%           - step                  : Convergence (relaxation) factor.
%           - filterOrderNo         : Order of the FIR filter.
%           - initialCoefficients   : Initial filter coefficients.  (COLUMN vector)
%           - gamma                 : Regularization factor.
%                                     (small positive constant to avoid singularity)
%           - wordlength            : wordlength for the fixed-point
%           - fractionLength        : number of bits after the binary point
%
%   Output Arguments:
%       . outputVector      :   Store the estimated output of each iteration.   (COLUMN vector)
%       . errorVector       :   Store the error for each iteration.             (COLUMN vector)
%       . coefficientVector :   Store the estimated coefficients for each iteration.
%                               (Coefficients at one iteration are COLUMN vector)
%
%   Authors:
%       . Guilherme de Oliveira Pinto   - guilhermepinto7@gmail.com & guilherme@lps.ufrj.br
%       . Markus Vinícius Santos Lima   - mvsl20@gmailcom           & markus@lps.ufrj.br
%       . Wallace Alves Martins         - wallace.wam@gmail.com     & wallace@lps.ufrj.br
%       . Luiz Wagner Pereira Biscainho - cpneqs@gmail.com          & wagner@lps.ufrj.br
%       . Paulo Sergio Ramirez Diniz    -                             diniz@lps.ufrj.br
%


%   Some Variables and Definitions:
%       . prefixedInput         :   Input is prefixed by nCoefficients -1 zeros.
%                                   (The prefix led to a more regular source code)
%
%       . regressor             :   Auxiliar variable. Store the piece of the
%                                   prefixedInput that will be multiplied by the
%                                   current set of coefficients.
%                                   (regressor is a COLUMN vector)
%
%       . nCoefficients         :   FIR filter number of coefficients.
%
%       . nIterations           :   Number of iterations.


%   Initialization Procedure
nCoefficients = S.filterOrderNo + 1;
nIterations = length(d);
word_len = S.wordlength;
frac_len = S.fractionLength;

%   Pre Allocations
e = fi(zeros(nIterations, 1), 1, word_len, frac_len);
y = fi(zeros(nIterations, 1), 1, word_len, frac_len);
w = fi(zeros(nCoefficients, (nIterations+1)), 1, word_len, frac_len);

%   Initial State Weight Vector
w(:,1) = fi(S.initCoeffs, 1, word_len, frac_len);

%   Improve source code regularity
prefixedInput = fi([zeros(nCoefficients-1,1)
                    transpose(x)], 1, word_len, frac_len);

%   Body
for it = 1:nIterations

    regressor =  prefixedInput(it+(nCoefficients-1):-1:it,1);

    y(it,1) = (w(:,it)')*regressor;

    e(it,1) = fi(d(it), 1, word_len, frac_len) - y(it,1);
    
    den = fi(S.gamma, 1, word_len, frac_len) + regressor'*regressor;
    prod = 1/den;
    num = regressor.*(fi(S.step, 1, word_len, frac_len)*e(it));
    
    w(:,it+1) = w(:,it) + num.*prod;

end

%   EOF
