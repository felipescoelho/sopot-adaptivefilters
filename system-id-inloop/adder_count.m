function [y, tc] = adder_count(h) 

%   adder_count.m
%       Counts the number of SOPOTs used in the implementation of the
%       filter given by 'h'.
%       This function acompains the mansucript "FIR Filter Design Based on 
%       Successive Approximation of Vectors" by Da Silva, Lovisolo, Dutra
%       and Diniz.
% 
%       Modified by Luiz Felipe Coelho in 23.11.2018, adusted for the usage
%       in SPT approximations for adaptive filters.
%
%       Remodified by Luiz Felipe Coelho in 17.05.2019, adjusted for the
%       usage in spt approximations for generic signals and numbers.
%
%   Syntax:
%       [y, tc] = adder_count(h)
% 
%
%   Output:
%       . y     : Total number of adders used for every filter component.
%       . tc    : Partial count for each filter component.
%
%
%   Input:
%       . h     : Filter whose SOPOTs are to be counted
%
%
%   Usage example:
%       Being h a FIR filter use
%
%           adder_count(h) 
%
%       To count how many SOPOTs are required in its implementation.


%   Observation: 
%       If the desired filter to count its adders is a symmetric or 
%       anti-symmetric filter, use just the rellevant (half) part.

%                   ##### Beginning the count: #####
count = 0 ;


%    ##### Makes sure the filter has a row-vector representation: #####
[nr, ~] = size(h);
if nr ~= 1
	h = h';
end
[~, nc] = size(h);  % nc: numbers to count


%           ##### Sets up array of possible powers-of-two: #####
h_pot = zeros(1, 129);


%     ##### Counts the number of terms for each vector component: #####
term_count = zeros(1, nc);
for i = 1:nc
	h_term = abs(h(i));
	residue = h_term;
	while residue ~= 0
		pot_term = round(log2(residue));
        h_pot(pot_term + 65) = 1;
		residue = abs(residue - 2^pot_term);
		term_count(i) = term_count(i) + 1;
		count = count + 1;
	end
end
y = count ;
tc = term_count ;
end



