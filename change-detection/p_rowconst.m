function z = p_rowconst(x, y, maxSumsfilt, maxPot)

%   p_rowconst.m
%       Implements the product of a row array with a constant where both
%       are represented in SOPOT. This function uses the function inner_p.m
%       ad it is necessary for the change detedetion.
%
%
%   Syntax:
%       z = p_rowconst(x, y, maxSumsfilt, maxPot)
%
%
%   Input args.:
%       . x     : array representing the row vector.
%                                     2-D ARRAY(bitplaneDepth, window_size)
%       . y     : array representing the constant.
%                                     ROW(1, bitplaneDepth)
%
%
%   Output args.:
%       . z     : array containing the resulting row.
%                                     2-D ARRAY(bitplaneDepth, window_size)
%
%
%   Author:
%       . Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
%
%

[word_len, win_size] = size(x);

% computing
z = zeros(word_len, win_size);
for i = 1:win_size
    z(:, i) = inner_p(transpose(x(:, i)), y, maxSumsfilt, maxPot);
end