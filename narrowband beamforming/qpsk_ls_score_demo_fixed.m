clc;
clear;
close all;
rng(7);

%% 1. Array geometry
M = 8;
lambda = 1;
d = lambda/2;
position = (0:M-1).' * d;

steering = @(theta) exp(-1j*2*pi/lambda * position * sind(theta));

theta_s = 20;                 % Desired signal DOA; used only to generate/evaluate
theta_i = -30;                % Interference DOA

a_s = steering(theta_s);
a_i = steering(theta_i);

%% 2. Desired pulse-shaped QPSK
Fs = 4e6;
Rs = 1e6;
sps = Fs/Rs;                  % 4 samples/symbol

rolloff = 0.35;
span = 10;
rrc = rcosdesign(rolloff, span, sps, 'sqrt');
gd = span*sps/2;

Nsym = 30000;                 % More samples reduce finite-sample cyclic leakage
guard = span + 5;

data_s_all = randi([0,3], Nsym + 2*guard, 1);
sym_s_all = pskmod(data_s_all, 4, pi/4, 'gray').';

s_full = upfirdn(sym_s_all, rrc, sps, 1);

Nsamp = Nsym*sps;
start_s = gd + guard*sps + 1;
s = s_full(start_s:start_s+Nsamp-1);

sym_s = sym_s_all(guard+1:guard+Nsym);

s = s / sqrt(mean(abs(s).^2));

%% 3. QPSK interference with a different symbol rate
sps_i = 5;
Ri = Fs/sps_i;                % 0.8 MHz

rrc_i = rcosdesign(rolloff, span, sps_i, 'sqrt');
gd_i = span*sps_i/2;
guard_i = span + 5;

Nsym_i = ceil(Nsamp/sps_i) + 2*guard_i + 5;
data_i = randi([0,3], Nsym_i, 1);
sym_i = pskmod(data_i, 4, pi/4, 'gray').';

i_full = upfirdn(sym_i, rrc_i, sps_i, 1);
start_i = gd_i + guard_i*sps_i + 1;
i_sig = i_full(start_i:start_i+Nsamp-1);
i_sig = i_sig / sqrt(mean(abs(i_sig).^2));

%% 4. Form array data
Ps = 1;
SNR_dB = 15;
INR_dB = 20;

Pn = Ps/10^(SNR_dB/10);
Pi = Pn*10^(INR_dB/10);

X_s = sqrt(Ps) * a_s * s;
X_i = sqrt(Pi) * a_i * i_sig;
V = sqrt(Pn/2) * (randn(M,Nsamp) + 1j*randn(M,Nsamp));

X = X_s + X_i + V;

%% 5. Select the desired cyclic frequency
% Desired signal: Rs = 1 MHz
% Interference:    Ri = 0.8 MHz
alpha = Rs;

% For this pulse shape, D = 0 gives a usable ordinary cyclic correlation.
% This is not a universal choice; it is chosen for this clean teaching case.
D = 0;

X_now = X(:,D+1:end);
X_delay = X(:,1:end-D);
K = size(X_now,2);
n = D:(Nsamp-1);

phase_pos = exp( 1j*2*pi*alpha/Fs*n);
phase_neg = exp(-1j*2*pi*alpha/Fs*n);

%% 6. Construct the LS-SCORE pseudo-reference
c = zeros(M,1);
c(1) = 1;
% c = rand(M,1);
% d_alpha[n] = c^H x[n-D] exp(j2*pi*alpha*n/Fs)
d_alpha = (c' * X_delay) .* phase_pos;

%% 7. Estimate Rxx and the input-reference correlation vector
Rxx = X_now*X_now'/K;

% p_xd = E{x[n] d_alpha^*[n]}
p_xd = X_now*d_alpha'/K;

% Equivalent cyclic-correlation form
Ralpha = (X_now .* phase_neg) * X_delay'/K;
p_xd_check = Ralpha*c;

fprintf('p_xd equivalence error = %.3e\n', ...
    norm(p_xd-p_xd_check)/max(norm(p_xd),eps));

%% 8. Diagnostic: does p_xd point toward the desired steering vector?
% This uses the true a_s only for simulation diagnosis, not for LS-SCORE.
alignment = abs(a_s' * p_xd)/(norm(a_s)*norm(p_xd));
fprintf('Alignment between p_xd and desired steering vector = %.4f\n', alignment);

%% 9. Solve the LS-SCORE Wiener equation
% delta = 1e-2 * real(trace(Rxx))/M;
% w_score = (Rxx + delta*eye(M)) \ p_xd;
w_score = Rxx \ p_xd;
% Blind power normalization; no desired DOA is used
% w_score = w_score / sqrt(real(w_score' * Rxx * w_score));

y_score = w_score' * X;

%% 10. Evaluate target/interference responses
% These quantities use known simulation parameters only for evaluation.
Gs = abs(w_score' * a_s)^2;
Gi = abs(w_score' * a_i)^2;
Gnoise = real(w_score' * w_score);

SINR_in = Ps/(Pi + Pn);
SINR_out = Ps*Gs/(Pi*Gi + Pn*Gnoise);

fprintf('Input SINR              = %.2f dB\n', 10*log10(SINR_in));
fprintf('LS-SCORE output SINR    = %.2f dB\n', 10*log10(SINR_out));
fprintf('Desired spatial gain    = %.4f\n', Gs);
fprintf('Interference spatial gain = %.4e\n', Gi);

%% 11. Beam pattern
theta_scan = -90:0.1:90;
A_scan = exp(-1j*2*pi/lambda * position * sind(theta_scan));
B = abs(w_score' * A_scan);

% Plot relative to the desired-direction response.
% This is only an evaluation normalization, so the desired direction is 0 dB.
B_target_ref_dB = 20*log10(B/(abs(w_score' * a_s) + eps) + 1e-12);

figure;
plot(theta_scan, B_target_ref_dB, 'LineWidth', 1.3);
hold on;
xline(theta_s, ':', 'Desired QPSK');
xline(theta_i, ':', 'QPSK interference');
grid on;
xlim([-90 90]);
ylim([-80 20]);
xlabel('DOA / degree');
ylabel('Response relative to desired direction / dB');
title('LS-SCORE beam pattern');
legend('LS-SCORE','Location','best');

%% 12. Matched filtering and constellation
x1_mf = conv(X(1,:), rrc);
y_mf = conv(y_score, rrc);

sample_index = gd + 1 + (0:Nsym-1)*sps;
x1_sym = x1_mf(sample_index);
y_sym = y_mf(sample_index);

valid = (span+1):(Nsym-span);
s_ref = sym_s(valid);
x1_valid = x1_sym(valid);
y_valid = y_sym(valid);

% Oracle complex-scalar alignment for plotting only
gain_x1 = (s_ref*x1_valid')/(x1_valid*x1_valid');
gain_y = (s_ref*y_valid')/(y_valid*y_valid');

x1_aligned = gain_x1*x1_valid;
y_aligned = gain_y*y_valid;

evm_x1 = sqrt(mean(abs(x1_aligned-s_ref).^2)/mean(abs(s_ref).^2));
evm_y  = sqrt(mean(abs(y_aligned-s_ref).^2)/mean(abs(s_ref).^2));

fprintf('First-sensor EVM        = %.2f %%\n', 100*evm_x1);
fprintf('LS-SCORE output EVM     = %.2f %%\n', 100*evm_y);

Nplot = 2500;
idx = max(1,length(valid)-Nplot+1):length(valid);

figure;
subplot(1,3,1);
plot(real(s_ref(idx)), imag(s_ref(idx)), '.');
axis equal; grid on;
xlabel('In-phase'); ylabel('Quadrature');
title('Transmitted QPSK symbols');

subplot(1,3,2);
plot(real(x1_aligned(idx)), imag(x1_aligned(idx)), '.');
axis equal; grid on;
xlabel('In-phase'); ylabel('Quadrature');
title('First sensor');

subplot(1,3,3);
plot(real(y_aligned(idx)), imag(y_aligned(idx)), '.');
axis equal; grid on;
xlabel('In-phase'); ylabel('Quadrature');
title('LS-SCORE output');

%% 13. Cyclic-frequency scan at the first sensor
alpha_grid = linspace(0.6e6,1.1e6,301);
cyclic_strength = zeros(size(alpha_grid));

x1_now = X(1,D+1:end);
x1_delay = X(1,1:end-D);

for k = 1:length(alpha_grid)
    ph = exp(-1j*2*pi*alpha_grid(k)/Fs*n);
    cyclic_strength(k) = abs(mean(x1_now .* conj(x1_delay) .* ph));
end

figure;
plot(alpha_grid/1e6, cyclic_strength, 'LineWidth',1.2);
hold on;
xline(Rs/1e6, ':', 'Desired 1.0 MHz');
xline(Ri/1e6, ':', 'Interference 0.8 MHz');
grid on;
xlabel('Cyclic frequency \alpha / MHz');
ylabel('|R_x^\alpha[D]|');
title('Cyclic-frequency selectivity');
