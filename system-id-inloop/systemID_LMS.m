% System Identification Problem
%   This script implements the SOPOT approximation in the LMS adaptive
%   filter.
% Paper:
%
% Authors:
%   Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
%   Lisandro Lovisolo - lisandro@uerj.br
%   Michel P. Tcheou - mtcheou@uerj.br
%

clc
clear
close all

% definitions
rng(42)
ensemble = 150;  % number of realizations
K = 10000;  % number of iterations
N = 4;  % number of coefficients
H = randn(N, 1);  % unknown system
H = H./norm(H);  % normalization
normH = norm(H);
mu_sopot = [.075 .075 .075 .075 .075 .075 .075 .075];  % step size (SOPOT)
mu_fi = .075;
max_sums = 8;  % maximum number of non-zero SPT terms
max_pot = 3;  % highest exponent for the SOPOT
word_len = 16;  % maximum number of SPT terms
word_len_fi = [3 5 7 9 11 13 15 16];  % wordlength for the fixed-point
frac_len_fi = [0 0 2 4 6 8 10 11];  % number of bits after the point, fixed-point
SNR = [-3 3 10 15];  % signal to noise ratio [dB]
snr_len = length(SNR);  % number of different SNR scenarios

% memory allocation
% fixed-point
w_fi = zeros(N, K+1, max_sums, snr_len, ensemble);  % coefficient vector (fixed-point)
norm_w_fi = zeros(K+1, max_sums, snr_len, ensemble);
MSE_fi = zeros(K, max_sums, snr_len, ensemble);
MSD_fi = zeros(K+1, max_sums, snr_len, ensemble);
% SOPOT
w_sopot = zeros(N, K+1, max_sums, snr_len, ensemble);  % coefficient vector (SOPOT)
norm_w_sopot = zeros(K+1, max_sums, snr_len, ensemble);
MSE_sopot = zeros(K, max_sums, snr_len, ensemble);
MSD_sopot = zeros(K+1, max_sums, snr_len, ensemble);

S = struct('filterOrderNo', N-1, 'initCoeffs', w_fi(:, 1, 1, 1, 1),...
    'wordLength', word_len, 'maxPot', max_pot);

% computing
for mc = 1:ensemble
    tic
    % some indication of the process
    fprintf('\n----------------------------------------\n')
    fprintf('Monte Carlo -- %u / %u \n', mc, ensemble)
    % creating signals
    x = randn(K, 1);  % input singal
    sigma_x2 = var(x);  % signal power
    x = x./sqrt(sigma_x2);  % signal power normalization
    d_aux = filter(H, 1, x);  % passing the signal through the "unknown" system
    Pd_aux = (norm(d_aux, 2).^2)./K;  % d_aux power
    n = randn(K, 1);  % some noise
    sigma_n2 = var(n);  % noise signal power
    for snr = 1:snr_len
        fprintf('SNR = %i dB\n', SNR(snr))
        Pn = Pd_aux*(10.^(-.1.*SNR(snr)));  % noise power, after adjustment
        n = sqrt(Pn./sigma_n2).*n;  % noise power adjustment
        d = d_aux + n;  % desired signal
        for spt = 1:max_sums
            fprintf('%u active SPT terms\n', spt)
            % setup for the fixed-point arithmetc
            G = globalfimath;
            G.ProductMode = 'SpecifyPrecision';
            G.SumMode = 'SpecifyPrecision';
            G = globalfimath('ProductWordLength', word_len_fi(spt),...
                'ProductFractionLength', frac_len_fi(spt),...
                'SumWordLength', word_len_fi(spt),...
                'SumFractionLength', frac_len_fi(spt));
            
            % fixed-point LMS
            S.step = mu_fi;
            S.wordlength = word_len_fi(spt);
            S.fractionLength = frac_len_fi(spt);
            [~, e_fi, w_fi(:, :, spt, snr, mc)] = LMS_fixed_point(d, x.', S);
            
            % SOPOT LMS
            S.step = mu_sopot(spt);
            S.maxSums = spt;
            S.maxSumsFilt = spt;
            [~, e_sopot, w_sopot(:, :, spt, snr, mc)] = LMS_spt(d, x, S);
            
            % MSE and MSD
            for k = 1:K+1
                norm_w_sopot(k, spt, snr, mc) = norm(w_sopot(:, k, spt, snr, mc));
                norm_w_fi(k, spt, snr, mc) = norm(single(w_fi(:, k, spt, snr, mc)));
            end
            deviation_sopot(:, 1) = normH - norm_w_sopot(:, spt, snr, mc);
            deviation_fi(:, 1) = normH - norm_w_fi(:, spt, snr, mc);
            MSD_fi(:, spt, snr, mc) = MSD_fi(:, spt, snr, mc) + abs(deviation_fi).^2;
            MSD_sopot(:, spt, snr, mc) = MSD_sopot(:, spt, snr, mc) + abs(deviation_sopot).^2;
            MSE_fi(:, spt, snr, mc) = MSE_fi(:, spt, snr, mc) + abs(e_fi(:, 1)).^2;
            MSE_sopot(:, spt, snr, mc) = MSE_sopot(:, spt, snr, mc) + abs(e_sopot(:, 1)).^2;
        end
    end
    if mc == 1
        elapsed = toc;
        sec = floor((elapsed./60 - floor(elapsed./60)).*60);
        min_aux = floor(elapsed./60);
        min = (min_aux./60 - floor(min_aux./60)).*60;
        hour_aux = floor(min_aux./60);
        hour = (hour_aux./24 - floor(hour_aux./24)).*24;
        day = floor(hour_aux./24);
        % time estimation
        estimated = ensemble.*elapsed;
        sec_est = floor((estimated./60 - floor(estimated./60)).*60);
        min_aux = floor(estimated./60);
        min_est = (min_aux./60 - floor(min_aux./60)).*60;
        hour_aux = floor(min_aux./60);
        hour_est = (hour_aux./24 - floor(hour_aux./24)).*24;
        day_est = floor(hour_aux./24);
    else
        elapsed = elapsed + toc;
        sec = floor((elapsed./60 - floor(elapsed./60)).*60);
        min_aux = floor(elapsed./60);
        min = (min_aux./60 - floor(min_aux./60)).*60;
        hour_aux = floor(min_aux./60);
        hour = (hour_aux./24 - floor(hour_aux./24)).*24;
        day = floor(hour_aux./24);
    end
    fprintf('The elapsed time is %u days, %u hours, %u minutes and %u seconds\n',...
        day, hour, min, sec);
    fprintf('The estimated time to finish is %u days, %u hours, %u minutes and %u seconds\n',...
        day_est, hour_est, min_est, sec_est);
end


% averaging
av_w_fi = sum(w_fi, 5)./ensemble;
av_w_sopot = sum(w_sopot, 5)./ensemble;
av_norm_w_fi = sum(norm_w_fi, 4)./ensemble;
av_norm_w_sopot = sum(norm_w_sopot, 4)./ensemble;
av_MSD_fi = sum(MSD_fi, 4)./ensemble;
av_MSD_sopot = sum(MSD_sopot, 4)./ensemble;
av_MSE_fi = sum(MSE_fi, 4)./ensemble;
av_MSE_sopot = sum(MSE_sopot, 4)./ensemble;

save('systemID_LMS.mat', 'av_*', 'H');

