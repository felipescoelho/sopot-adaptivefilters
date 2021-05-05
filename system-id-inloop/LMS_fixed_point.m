function    [y, e, w] =   LMS_fixed_point(d, x, S)

%   LMS_fixed_point.m
%       Implements the Complex LMS algorithm for COMPLEX valued data.
%       (Algorithm 3.2 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
%
% - Update: 25.05.2020 (dd.mm.yyyy)
%   Modified by Luiz Felipe da S. Coelho to be used as a fixed-point LMS
%                                    - email: luizfelipe.coelho@smt.ufrj.br
%
%   Syntax:
%       [y, e, w] = LMS_fixed_point(d, x, S)
%
%   Input Arguments:
%       . d         : Desired signal.                               (ROW vector)
%       . x         : Signal fed into the adaptive filter.          (ROW vector)
%       . S         : Structure with the following fields
%           - step              : Convergence (relaxation) factor.
%           - filterOrderNo     : Order of the FIR filter.
%           - initCoeffs        : Initial filter coefficients.  (COLUMN vector)
%           - wordlength        : Wordlength for the fixed-point
%           - fractionlength    : number of bits after the binary point
%
%   Output Arguments:
%       . y : Store the estimated output of each iteration.   (COLUMN vector)
%       . e : Store the error for each iteration.             (COLUMN vector)
%       . w : Store the estimated coefficients for each iteration.
%                               (Coefficients at one iteration are COLUMN vector)
%
%   Authors:
%       . Guilherme de Oliveira Pinto   - guilhermepinto7@gmail.com & guilherme@lps.ufrj.br
%       . Markus Vinícius Santos Lima   - mvsl20@gmailcom           & markus@lps.ufrj.br
%       . Wallace Alves Martins         - wallace.wam@gmail.com     & wallace@lps.ufrj.br
%       . Luiz Wagner Pereira Biscainho - cpneqs@gmail.com          & wagner@lps.ufrj.br
%       . Paulo Sergio Ramirez Diniz    -                             diniz@lps.ufrj.br
%
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
nCoefficients       =   S.filterOrderNo+1;
nIterations         =   length(d);
word_len = S.wordlength;
frac_len = S.fractionLength;

%   Pre Allocations
e = fi(zeros(nIterations   ,1), 1, word_len, frac_len);
y = fi(zeros(nIterations   ,1), 1, word_len, frac_len);
w = fi(zeros(nCoefficients ,(nIterations+1)), 1, word_len, frac_len);

%   Initial State Weight Vector
w(:,1)  =   fi(S.initCoeffs, 1, word_len, frac_len);

%   Improve source code regularity
prefixedInput           =   fi([zeros(nCoefficients-1,1)
                             transpose(x)], 1, word_len, frac_len);

%   Body
for it = 1:nIterations

    regressor = prefixedInput(it+(nCoefficients-1):-1:it,1);

    y(it,1)          =   (w(:,it)')*regressor;

    e(it,1)           =   fi(d(it), 1, word_len, frac_len)-y(it,1);

    w(:,it+1) = w(:,it)+(fi(S.step, 1, word_len, frac_len)*...
        conj(e(it,1))*regressor);

end

%   EOF
