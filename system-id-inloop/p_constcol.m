function z = p_constcol(x, y, maxs_filt, max_pot)

% p_constcol.m
%

[N, word_len] = size(y);

z = zeros(N , word_len);

for n = 1:N
    z(n, :) = inner_p(x, y(n, :), maxs_filt, max_pot);
end