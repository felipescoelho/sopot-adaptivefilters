% System Identification Problem, using RLS_spt
% Used for the results in the ISCAS 2020 paper.
% author: Luiz Felipe da S. Coelho -- lfscoelho@ieee.org

clear
close all
clc


% Definitions:
rng(5555);
ensemble = 1000;  % number of realizations
K = 10000;  % number of iterations
H = randn(4, 1);  % unknown system
normH = norm(H);  % unknown system's norm
H = H/normH;  % unknown system normalization
sigma_n2 = 1e-2;  % noise power
sigma_n2_2 = .5;
N = 4;  % number of coefficients
delta = .09;  % small positive constant (to initialize S_d)
lambda = [.975 .975 .96 .96 .955 .9535 .9534 .954];  % forgetting factor
lambda2 = .9987;
maxSums = 8;  % maximum number of non-zero terms for the spt approximtion
maxPot = 3;  % the greater power for the spt (1, 2 or 3)
bitplaneDepth = 16;  % maximum bit-plane depth
wordlength = [3 5 7 9 11 13 15 16];  % wordlength for the fixed-point, in bits
frac_len = [0 0 2 4 6 8 10 11];  % Depth of the numeric representation


%   Allocating memory:
W = zeros(N, (K+1), maxSums, ensemble);  % Coefficient vector (each iteration)
W_spt = zeros(N, (K+1), maxSums, ensemble);  % coefficients for the spt filter
MSE = zeros(K, maxSums, ensemble);  % MSE
MSE_spt = zeros(K, maxSums, ensemble);  % MSE for the spt filtering
MSEmin = zeros(K, maxSums, ensemble);  % minimal MSE (each iteration)
MSD = zeros(K+1, maxSums, ensemble);
MSD_spt = zeros(K+1, maxSums, ensemble);
W_norm_aux = zeros(K+1, maxSums, ensemble);
W_spt_norm_aux = zeros(K+1, maxSums, ensemble);
W_aux = zeros(K+1, maxSums, ensemble);
W_spt_aux = zeros(K+1, maxSums, ensemble);
e_depth = zeros(K, maxSums);
y_depth = zeros(K, maxSums);
desired_depth = zeros(K, maxSums);
w_depth = zeros(K, N, maxSums);


%   Computing:
fprintf('System Identification using the SPT-RLS:\n')
fprintf(' SNR = 10 dB \n')
tic;
for l = 1:ensemble
    fprintf('\n--------------------------------------------\n')
    fprintf('Monte Carlo number: %u \n\n', l)
    X = zeros(N, 1);  % input in certain iteration
    d = zeros(1, K);  % desired signal
    x = randn(K, 1);  % input noise
    sigma_x2 = var(x);  % input noise power
    n = sqrt(sigma_n2)*randn(K,1);  % a noise
        
    for k = 1: K
        X = [x(k,1)
            X(1:(N-1), 1)];  % input signal (tapped delay line)
        d(k) = (H'*X(:, 1)) + n(k);  % desired signal
    end
    
    for spt_terms = 1:maxSums
        fprintf('\n------------------------------\n')
        fprintf('For %u SPT terms: \n\n', spt_terms)
        % setup for the fixed-point arithmetic
        G = globalfimath;
        G.ProductMode = 'SpecifyPrecision';
        G.SumMode = 'SpecifyPrecision';
        G = globalfimath('ProductWordLength', wordlength(spt_terms),...
            'ProductFractionLength', frac_len(spt_terms),...
            'SumWordLength', wordlength(spt_terms),...
            'SumFractionLength', frac_len(spt_terms));
        
        % Fixed-point adaptive filter
        S = struct('filterOrderNo', (N-1), 'initialCoefficients',...
            W(:, 1, l), 'delta', delta, 'lambda', lambda2,...
            'wordlength', wordlength(spt_terms), 'fractionLength',...
            frac_len(spt_terms));
        [y, e, W(:, :, spt_terms, l)] = RLS_Alt_fixed_point_2(d, x.', S);
        
        % SPT-RLS adaptive filter
        S_spt = struct('filterOrder', (N-1), 'initCoeff', W(:, 1, 1),...
            'delta', delta, 'lambda', lambda(spt_terms),...
            'bitplaneDepth', bitplaneDepth, 'maxSums', spt_terms,...
            'maxSumsfilt', spt_terms, 'maxPot', maxPot);
        [y_spt, e_spt, W_spt(:, :, spt_terms, l),...
            y_depth(:, spt_terms), e_depth(:, spt_terms),...
            desired_depth(:, spt_terms), w_depth(:, :, spt_terms)] = RLS_spt(d.', x, S_spt);
        
        %   MSE/MSD calculations
        for iter = 1:K+1
            W_aux(iter, spt_terms, l) = norm(H - W(:, iter, spt_terms, l));
            W_spt_aux(iter, spt_terms, l) = norm(H - W_spt(:, iter, spt_terms, l));
        end
        for iter = 1:K+1
            W_norm_aux(iter, spt_terms, l) = norm(single(W(:, iter, spt_terms, l)));
            W_spt_norm_aux(iter, spt_terms, l) = norm(W_spt(:, iter, spt_terms, l));
        end
        MSD(:, spt_terms, l) = MSD(:, spt_terms, l) + (abs(W_aux(:, spt_terms, l))).^2;
        MSD_spt(:, spt_terms, l) = MSD_spt(:, spt_terms, l) + (abs(W_spt_aux(:, spt_terms, l))).^2;
        MSE_spt(:, spt_terms, l) = MSE_spt(:, spt_terms, l) + (abs(e_spt(:, 1))).^2;
        MSE(:, spt_terms, l) = MSE(:, spt_terms, l) + (abs(e(:, 1))).^2;
        MSEmin(:, spt_terms, l) = MSEmin(:, spt_terms, l)+(abs(n(:))).^2;
    end
end
toc()


%  SNR = 3

rng(5555);
W_2 = zeros(N, (K+1), maxSums, ensemble);  % Coefficient vector (each iteration)
W_spt_2 = zeros(N, (K+1), maxSums, ensemble);  % coefficients for the spt filter
MSE_2 = zeros(K, maxSums, ensemble);  % MSE
MSE_spt_2 = zeros(K, maxSums, ensemble);  % MSE for the spt filtering
MSEmin_2 = zeros(K, maxSums, ensemble);  % minimal MSE (each iteration)
MSD_2 = zeros(K+1, maxSums, ensemble);
MSD_spt_2 = zeros(K+1, maxSums, ensemble);
W_norm_aux_2 = zeros(K+1, maxSums, ensemble);
W_spt_norm_aux_2 = zeros(K+1, maxSums, ensemble);
W_aux_2 = zeros(K+1, maxSums, ensemble);
W_spt_aux_2 = zeros(K+1, maxSums, ensemble);

fprintf(' SNR = 3 dB \n')
for l = 1:ensemble
    fprintf('\n--------------------------------------------\n')
    fprintf('Monte Carlo repetition number: %u \n\n', l)
    X = zeros(N, 1);  % input in certain iteration
    d = zeros(1, K);  % desired signal
    x = randn(K, 1);  % input noise
    sigma_x2 = var(x);  % input noise power
    n = sqrt(sigma_n2_2)*randn(K,1);  % a noise
        
    for k = 1: K
        X = [x(k,1)
            X(1:(N-1), 1)];  % input signal (tapped delay line)
        d(k) = (H'*X(:, 1)) + n(k);  % desired signal
    end
    
    for spt_terms = 1:maxSums
        fprintf('\n------------------------------\n')
        fprintf('For %u SPT terms: \n\n', spt_terms)
        % setup for the fixed-point arithmetic
        G = globalfimath;
        G.ProductMode = 'SpecifyPrecision';
        G.SumMode = 'SpecifyPrecision';
        G = globalfimath('ProductWordLength', wordlength(spt_terms),...
            'ProductFractionLength', frac_len(spt_terms),...
            'SumWordLength', wordlength(spt_terms),...
            'SumFractionLength', frac_len(spt_terms));
        
        % Fixed-point adaptive filter
        S = struct('filterOrderNo', (N-1), 'initialCoefficients',...
            W_2(:, 1, l), 'delta', delta, 'lambda', lambda2,...
            'wordlength', wordlength(spt_terms), 'fractionLength',...
            frac_len(spt_terms));
        [y, e, W_2(:, :, spt_terms, l)] = RLS_Alt_fixed_point_2(d, x.', S);
        
        % SPT-RLS adaptive filter
        S_spt = struct('filterOrder', (N-1), 'initCoeff', W_2(:, 1, 1),...
            'delta', delta, 'lambda', lambda(spt_terms),...
            'bitplaneDepth', bitplaneDepth, 'maxSums', spt_terms,...
            'maxSumsfilt', spt_terms, 'maxPot', maxPot);
        [y_spt, e_spt, W_spt_2(:, :, spt_terms, l),...
            y_depth(:, spt_terms), e_depth(:, spt_terms),...
            desired_depth(:, spt_terms), w_depth(:, :, spt_terms)] = RLS_spt(d.', x, S_spt);
        
        %   MSE/MSD calculations
        for iter = 1:K+1
            W_aux(iter, spt_terms, l) = norm(H - W_2(:, iter, spt_terms, l));
            W_spt_aux(iter, spt_terms, l) = norm(H - W_spt_2(:, iter, spt_terms, l));
        end
        for iter = 1:K+1
            W_norm_aux_2(iter, spt_terms, l) = norm(single(W_2(:, iter, spt_terms, l)));
            W_spt_norm_aux_2(iter, spt_terms, l) = norm(W_spt_2(:, iter, spt_terms, l));
        end
        MSD_2(:, spt_terms, l) = MSD_2(:, spt_terms, l) + (abs(W_aux(:, spt_terms, l))).^2;
        MSD_spt_2(:, spt_terms, l) = MSD_spt_2(:, spt_terms, l) + (abs(W_spt_aux(:, spt_terms, l))).^2;
        MSE_spt_2(:, spt_terms, l) = MSE_spt_2(:, spt_terms, l) + (abs(e_spt(:, 1))).^2;
        MSE_2(:, spt_terms, l) = MSE_2(:, spt_terms, l) + (abs(e(:, 1))).^2;
        MSEmin_2(:, spt_terms, l) = MSEmin_2(:, spt_terms, l)+(abs(n(:))).^2;
    end
end

%   Averaging:
sav_W_av = sum(W, 4)/ensemble;
sav_W_spt_av = sum(W_spt, 4)/ensemble;
sav_w_norm = sum(W_norm_aux, 3)./ensemble;
sav_w_spt_norm = sum(W_spt_norm_aux, 3)./ensemble;
sav_MSD_av = sum(MSD, 3)./ensemble;
sav_MSD_spt_av = sum(MSD_spt, 3)./ensemble;
sav_MSE_av = sum(MSE, 3)/ensemble;
sav_MSE_spt_av = sum(MSE_spt, 3)/ensemble;
sav_MSEmin_av = sum(MSEmin, 3)/ensemble;

% SNR = 3
sav_W_av_2 = sum(W_2, 4)/ensemble;
sav_W_spt_av_2 = sum(W_spt_2, 4)/ensemble;
sav_w_norm_2 = sum(W_norm_aux_2, 3)./ensemble;
sav_w_spt_norm_2 = sum(W_spt_norm_aux_2, 3)./ensemble;
sav_MSD_av_2 = sum(MSD_2, 3)./ensemble;
sav_MSD_spt_av_2 = sum(MSD_spt_2, 3)./ensemble;
sav_MSE_av_2 = sum(MSE_2, 3)/ensemble;
sav_MSE_spt_av_2 = sum(MSE_spt_2, 3)/ensemble;
sav_MSEmin_av_2 = sum(MSEmin_2, 3)/ensemble;


%   Coefficient Evaluation:
% Coefficient Norm
sav_normH = norm(H)*ones(K, 1);

save('SystemID.mat','sav_*');

% figure(1),
% subplot(211),
% plot(1:300, normH(1:300), '--k'), hold on,
% plot(1:300, w_spt_norm(1:300, end-5)),
% plot(1:300, w_spt_norm(1:300, end-3)),
% plot(1:300, w_spt_norm(1:300, end)),
% plot(1:300, w_norm(1:300, end), ':'), grid on, hold off,
% title('(a)'),
% subplot(212),
% plot(1:300, normH(1:300), '--k'), hold on,
% plot(1:300, w_spt_norm_2(1:300, end-5)),
% plot(1:300, w_spt_norm_2(1:300, end-3)),
% plot(1:300, w_spt_norm_2(1:300, end)),
% plot(1:300, w_norm_2(1:300, end), ':'), grid on, hold off,
% title('(b)'),
% legend(' H', ' M/N = 3', ' M/N = 5', ' M/N = 8', ' 16-bit fixed-point'),
% xlabel('Iterations'), ylabel('Coeficients Norm')
% 
% % Absolute Error:
% figure(2),
% subplot(211),
% ax = gca;
% plot(1:300, 10*log10(MSD_spt_av(1:300, end-5))), hold on,
% plot(1:300, 10*log10(MSD_spt_av(1:300, end-4))),
% plot(1:300, 10*log10(MSD_spt_av(1:300, end-3))),
% ax.ColorOrderIndex = 1;
% plot(1:300, 10*log10(MSD_av(1:300, end-5)), ':'),
% plot(1:300, 10*log10(MSD_av(1:300, end-4)), ':'),
% plot(1:300, 10*log10(MSD_av(1:300, end-3)), ':'), grid on, hold off,
% title('(a)'),
% subplot(212),
% ax = gca;
% plot(1:300, 10*log10(MSD_spt_av_2(1:300, end-5))), hold on,
% plot(1:300, 10*log10(MSD_spt_av_2(1:300, end-4))),
% plot(1:300, 10*log10(MSD_spt_av_2(1:300, end-3))),
% ax.ColorOrderIndex = 1;
% plot(1:300, 10*log10(MSD_av_2(1:300, end-5)), ':'),
% plot(1:300, 10*log10(MSD_av_2(1:300, end-4)), ':'),
% plot(1:300, 10*log10(MSD_av_2(1:300, end-3)), ':'), grid on, hold off,
% title('(b)'),
% legend(' M/N = 3', ' M/N = 4', ' M/N = 5', ' 7-bit fixed-point', ' 9-bit fixed-point', ' 11-bit fixed-point'),
% xlabel('Iterations'), ylabel('MSD [dB]')
% 
% figure(3),
% subplot(211),
% plot(1:300, 10*log10(MSD_spt_av(1:300, end-5))), hold on,
% plot(1:300, 10*log10(MSD_spt_av(1:300, end-3))),
% plot(1:300, 10*log10(MSD_spt_av(1:300, end))),
% plot(1:300, 10*log10(MSD_av(1:300, end)), ':'), grid on, hold off,
% title('(a)'),
% subplot(212),
% plot(1:300, 10*log10(MSD_spt_av_2(1:300, end-5))), hold on,
% plot(1:300, 10*log10(MSD_spt_av_2(1:300, end-3))),
% plot(1:300, 10*log10(MSD_spt_av_2(1:300, end))),
% plot(1:300, 10*log10(MSD_av_2(1:300, end)), ':'), grid on, hold off,
% title('(b)'),
% legend(' M/N = 3', ' M/N = 5', ' M/N = 8', ' 16-bit fixed-point'),
% xlabel('Iterations'), ylabel('MSD [dB]')

%   MSE Plotting:
%{
figure(3),
subplot(311),
plot(1:K, 10*log10(MSE_av(:, end))), grid on
title('(a)')
subplot(312),
plot(1:K, 10*log10(MSE_spt_av(:, end))), grid on
title('(b)'), ylabel('MSE[dB]')
subplot(313),
plot(1:K, 10*log10(MSEmin_av(:, end))), grid on
title('(c)'), xlabel('Iterations'), ylabel('MSE[dB]')
%}
% figure(4),
% subplot(121),
% surf(10*log10(MSE_spt_av(1:300, :, :))), colormap(spring)
% title('(a)'),
% xlabel('M/N'), ylabel('Iterations'), zlabel('MSE[dB]')
% subplot(122),
% surf(10*log10(MSE_spt_av_2(1:300, :, :))), colormap(spring)
% title('(b)'),
% xlabel('M/N'), ylabel('Iterations'), zlabel('MSE[dB]')

%{
%   Plotting:
figure,
plot(1:K,10*log10(MSE_av),'-k');
title('Learning Curve for MSE');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

figure,
plot(1:K, 10*log10(MSE_spt_av), '-k');
title('Learning Curve for MSE_{SPT}');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');


figure,
plot(1:K,10*log10(MSEmin_av),'-k');
title('Learning Curve for MSEmin');
xlabel('Number of iterations, k'); ylabel('MSEmin [dB]');
%}


%   EoF


