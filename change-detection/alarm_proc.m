function [CD_time, fp_ratio, fn_ratio] = alarm_proc(alarm, ref, trans_time)

% alarm_proc.m
%
%   Function to identify the time for change detection (CD), and the false
%   positivve and negative ratio. Between the alarm signal and a referrence
%   signal, with a give transition time.
%
%   Input:
%       . alarm     : alarm binary mask (empty = true, occupied = false)
%       . ref       : referrence binary mask (empty = true, occupied = false)
%       . trans_time: accepted delay time for change detection
%
%   Output:
%       . CD_time   : time of CD  ROW VECTOR
%       . fp_ratio  : false positive ratio
%       . fn_ratio  : false negative ratio
%
%   Author:
%       . Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
%
%   Obs.:
%       We test for empty channel, i.e., when the channel is empty it
%       results in a positive.

sig_len = length(alarm);  % signal length (samples)
nBlocks = 9;  % we have 9 blocks of silence/signal, odd beeing silent
block_len = floor(sig_len./nBlocks);  % block length  (samples)
nSubBlocks = floor(block_len./trans_time);
pos_time = nSubBlocks*5;  % amount of silent samples
neg_time = nSubBlocks*4;  % amount of signal samples

CD_time = zeros(1, nBlocks-1);  % transition time counter
alarm_aux = ones(1, trans_time);
false_pos = 0;  % false positive count
false_neg = 0;  % false negative count


for block = 1:nBlocks
    first_blood = 0;  % token for first detection in a block
    sig_block = alarm((((block-1)*block_len)+1):(block*block_len));
    ref_block = ref((((block-1)*block_len)+1):(block*block_len));
    if block ~= 1  % detection of first alarm in a block
        for i = 1:80
            if sig_block(i) == ref_block(i) && first_blood == 0
                first_blood = 1;  % token for first detection
                CD_time(block-1) = i;  % time of CD in samples
            end
        end
    end
    if rem(block, 2) ~= 0  % odd block
        % This should be a silent block (all should be 1's)
        for subBlock = 1:nSubBlocks
            sig_subblock = sig_block((((subBlock-1)*trans_time)+...
                1):(subBlock*trans_time));  % this should be filled with true
            test = alarm_aux(sig_subblock);
            if isempty(test)
                % if is empty sig_subblock does not have any true
                false_neg = false_neg + 1;
            end
        end
    else  % even block
        % This should be a signal block, therefore all false
        for subBlock = 1:nSubBlocks
            sig_subblock = sig_block((((subBlock-1)*trans_time)+...
                1):(subBlock*trans_time));  % this should be filled with false
            test = alarm_aux(sig_subblock);
            if ~isempty(test)
                % if is not empty sig_subblock has a true
                false_pos = false_pos + 1;
            end
        end
    end
end

% false positive ratio = false positive / total negative
% false negative ratio = false negative / total positive

fp_ratio = false_pos./neg_time;
fn_ratio = false_neg./pos_time;







