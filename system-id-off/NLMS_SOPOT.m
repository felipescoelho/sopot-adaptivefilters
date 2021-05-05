function    [sopot_outVector,...
             sopot_errorVector,...
             sopotVector,...
             km] =   NLMS_SOPOTinloop(desired,input,S)

%   NLMS.m
%       Implements the Normalized LMS algorithm for COMPLEX valued data.
%       (Algorithm 4.3 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
%
%   Syntax:
%       [outputVector,errorVector,coefficientVector,sopotOutputVector, sopotErrorVector, sopotCoefficientVector] = NLMS(desired,input,S)
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
%           - adders                : Number of SOPOTs in the mpgbp
%                                     approximation
%
%   Output Arguments:
%       . outputVector      :   Store the estimated output of each iteration.   (COLUMN vector)
%       . errorVector       :   Store the error for each iteration.             (COLUMN vector)
%       . coefficientVector :   Store the estimated coefficients for each iteration.
%                               (Coefficients at one iteration are COLUMN vector)
%       . sopot_outVector   :   Store the estimated output for the SOPOT
%                               approximation of each iteration.                (COLUMN vector)
%       . sopot_errorVector :   Store the error for the SOPOT approximation
%                               for each iteration.                             (COLUMN vector)
%       . sopotVector       :   Store the estimated coefficients for the
%                               SOPOT approximation after each iteration.       (Coefficients at one iteration are COLUMN vector)
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
nCoefficients       =   S.filterOrderNo+1;
nIterations         =   length(desired);
nAdders             =   S.adders;
P                   =   floor(sqrt(S.filterOrderNo+1));

%   Pre Allocations
coefficientVector       =   zeros(nCoefficients ,(nIterations+1));
sopot_outVector         =   zeros(nIterations   ,1);
sopot_errorVector       =   zeros(nIterations   ,1);
sopotVector             =   zeros(nCoefficients ,(nIterations+1));

%   Initial State Weight Vector
coefficientVector(:,1)  =   S.initialCoefficients;
sopotVector(:,1)        =   S.initialCoefficients;

%   Improve source code regularity
prefixedInput           =   [zeros(nCoefficients-1,1)
                             transpose(input)];

%   Body
for it = 1:nIterations

    regressor                   =   prefixedInput(it+(nCoefficients-1):-1:it,1);

    sopot_outVector(it, 1)      =   (sopotVector(:,it)')*regressor;
    
    sopot_errorVector(it,1)     =   desired(it)-sopot_outVector(it,1);
    
    coefficientVector(:,it+1)   =   sopotVector(:,it)+(...
                                    (S.step/(S.gamma+regressor'*regressor))*...
                                    conj(sopot_errorVector(it,1))*regressor);

%   SOPOT-Cycle                                
    
    [h_spt, km]                 =   mpgbp_filter_approxP(coefficientVector(:,it+1), nAdders, P); 

    sopotVector(:,it+1)         =   (h_spt)';
   
                                
end

%   EOF
