clc;
clear;
close all;
rng(7);

%% Worst-case robust adaptive beamforming using CVX
% Solve:
%   min_w  w^H Rhat w
%   s.t.   Re{w^H a0} >= 1 + epsilon*||w||_2
%          Im{w^H a0} = 0
%
% Run cvx_setup before executing this script.

%% 1. Array and signal parameters
M = 12;                         % Number of array elements
N = 5e3;                        % Number of snapshots
m = (0:M-1).';                  % Element indices
d_over_lambda = 0.5;            % Element spacing / wavelength

theta0 = 0;                     % Presumed target direction (deg)
theta_s_true = 2;               % True target direction (deg)
theta_i = [-30, 40];            % Interference directions (deg)

SNR_dB = 10;
INR_dB = [30, 25];
sigma2 = 1;

Ps = sigma2 * 10^(SNR_dB/10);
Pi = sigma2 * 10.^(INR_dB/10);

% ULA steering vector; angles are measured from broadside.
steer = @(theta) exp(-1j*2*pi*d_over_lambda*m*sind(theta));

a0 = steer(theta0);
a_true = steer(theta_s_true);
a_i1 = steer(theta_i(1));
a_i2 = steer(theta_i(2));

%% 2. Convert the angular interval to a steering-vector error radius
Delta_theta = 2.5;              % Uncertainty half-width (deg)
theta_unc = linspace(theta0-Delta_theta, theta0+Delta_theta, 2001);
A_unc = steer(theta_unc);
epsilon = max(vecnorm(A_unc-a0, 2, 1));

fprintf('||a0||_2 = %.6f\n', norm(a0));
fprintf('epsilon  = %.6f\n', epsilon);
fprintf('||a_true-a0||_2 = %.6f\n\n', norm(a_true-a0));

if epsilon >= norm(a0)
    error('epsilon must be smaller than ||a0||; otherwise the robust constraint is infeasible.');
end

%% 3. Generate received data
% QPSK target signal
data = randi([0, 3], 1, N);
s_unit = exp(1j*(pi/4 + pi/2*data));
s = sqrt(Ps)*s_unit;

% Independent complex Gaussian interference signals
i1 = sqrt(Pi(1)/2)*(randn(1, N) + 1j*randn(1, N));
i2 = sqrt(Pi(2)/2)*(randn(1, N) + 1j*randn(1, N));

% Spatially white complex Gaussian noise
noise = sqrt(sigma2/2)*(randn(M, N) + 1j*randn(M, N));

% Received array data
X = a_true*s + a_i1*i1 + a_i2*i2 + noise;

% Raw sample covariance matrix. Ideal MVDR uses this matrix directly,
% without Hermitian smoothing or diagonal loading.
Rhat_raw = X*X'/N;

% Numerically stabilized covariance used by mismatched MVDR and WC-RAB.
Rhat = Rhat_raw;
% Rhat = (Rhat + Rhat')/2;

% Small diagonal loading for numerical stability
% delta_num = 1e-8*real(trace(Rhat))/M;
% Rhat = Rhat + delta_num*eye(M);

%% 4. Three beamformers
% 4.1 Ideal MVDR: use the true steering vector and the raw covariance.
q = Rhat_raw\a_true;
w_ideal = q/(a_true'*q);

% 4.2 Mismatched MVDR: use the presumed steering vector.
q = Rhat\a0;
w_mvdr = q/(a0'*q);

% 4.3 Worst-case robust adaptive beamformer
if exist('cvx_begin', 'file') ~= 2
    error('CVX is not available on the MATLAB path. Install CVX and run cvx_setup first.');
end

% Normalize the covariance before factorization. Multiplying the objective
% by a positive scalar does not change its minimizer, and this improves
% the numerical conditioning of the conic solver.

% Rhat_scale = real(trace(Rhat))/M;
% C = chol(Rhat/Rhat_scale);
C = chol(Rhat);
% Rhat/Rhat_scale = C'*C, so the objective is ||C*w||_2^2.
% This form also avoids a vec/quad_form compatibility issue in some
% CVX 2.2 installations used with recent MATLAB releases.
% vec =1 
cvx_clear
cvx_begin quiet
    variable w_wc(M) complex

    % minimize(w_wc' * Rhat * w_wc)
    minimize(sum_square_abs(C*w_wc))
    subject to
        real(w_wc' * a0) >= 1 + epsilon * norm(w_wc, 2);
        imag(w_wc' * a0) == 0;
cvx_end

fprintf('CVX status = %s\n', cvx_status);
fprintf('CVX optimum value = %.6e\n\n', cvx_optval);

if isempty(strfind(cvx_status, 'Solved')) %#ok<STREMP>
    error('CVX did not solve the robust beamforming problem. Status: %s', cvx_status);
end

%% 5. Verify the robust constraint
nominal_gain = w_wc'*a0;
wc_lower_bound = real(nominal_gain) - epsilon*norm(w_wc);

fprintf('Re{w_wc^H a0} = %.9f\n', real(nominal_gain));
fprintf('Im{w_wc^H a0} = %.3e\n', imag(nominal_gain));
fprintf('Re{w_wc^H a0}-epsilon||w_wc|| = %.9f\n\n', wc_lower_bound);

%% 6. Theoretical output SINR
Rin = Pi(1)*(a_i1*a_i1') + ...
      Pi(2)*(a_i2*a_i2') + ...
      sigma2*eye(M);

output_sinr = @(w, a) ...
    10*log10(Ps*abs(w'*a)^2 / real(w'*Rin*w));

SINR_ideal = output_sinr(w_ideal, a_true);
SINR_mvdr  = output_sinr(w_mvdr, a_true);
SINR_wc    = output_sinr(w_wc, a_true);

fprintf('Ideal MVDR output SINR      = %8.3f dB\n', SINR_ideal);
fprintf('Mismatched MVDR output SINR = %8.3f dB\n', SINR_mvdr);
fprintf('WC-RAB output SINR          = %8.3f dB\n', SINR_wc);
fprintf('|w_mvdr^H a_true| = %.6f\n', abs(w_mvdr'*a_true));
fprintf('|w_wc^H a_true|   = %.6f\n\n', abs(w_wc'*a_true));

%% 7. Beampattern
theta_scan = -90:0.1:90;
A_scan = steer(theta_scan);

B_ideal = abs(w_ideal'*A_scan);
B_mvdr  = abs(w_mvdr'*A_scan);
B_wc    = abs(w_wc'*A_scan);

B_ideal_dB = 20*log10(B_ideal/max(B_ideal) + 1e-12);
B_mvdr_dB  = 20*log10(B_mvdr/max(B_mvdr) + 1e-12);
B_wc_dB    = 20*log10(B_wc/max(B_wc) + 1e-12);

figure;
plot(theta_scan, B_ideal_dB, '--', 'LineWidth', 1);
hold on;
plot(theta_scan, B_mvdr_dB, '-.', 'LineWidth', 1);
plot(theta_scan, B_wc_dB, '-', 'LineWidth', 1);

xline(theta0, ':', 'Presumed SOI','LineWidth',1.5);
xline(theta_s_true, ':', 'True SOI','LineWidth',1.5);
xline(theta_i(1), ':', 'Interference 1','LineWidth',1.5);
xline(theta_i(2), ':', 'Interference 2','LineWidth',1.5);

grid on;
xlim([-90, 90]);
ylim([-80, 5]);
xlabel('Direction / degree');
ylabel('Normalized response / dB');
title('Worst-case robust adaptive beamforming');
legend('Ideal MVDR', 'Mismatched MVDR', 'WC-RAB', ...
       'Location', 'southoutside', 'NumColumns', 3);

%% 8. Scan the true target direction and compare robustness
theta_error = -5:0.25:5;

SINR_ideal_curve = zeros(size(theta_error));
SINR_mvdr_curve  = zeros(size(theta_error));
SINR_wc_curve    = nan(size(theta_error));

for k = 1:numel(theta_error)
    theta_k = theta0 + theta_error(k);
    a_k = steer(theta_k);

    Xk = a_k*s + a_i1*i1 + a_i2*i2 + noise;
    % Ideal MVDR uses the raw sample covariance without smoothing/loading.
    Rk_raw = Xk*Xk'/N;

    % Stabilized covariance for mismatched MVDR and WC-RAB.
    Rk = Rk_raw;
    Rk = (Rk + Rk')/2;
    Rk = Rk + 1e-8*real(trace(Rk))/M*eye(M);
    Rk_scale = real(trace(Rk))/M;
    Ck = chol(Rk/Rk_scale);

    % Ideal MVDR
    q = Rk_raw\a_k;
    w_ideal_k = q/(a_k'*q);

    % Mismatched MVDR
    q = Rk\a0;
    w_mvdr_k = q/(a0'*q);

    % WC-RAB
    cvx_begin quiet
        variable w_wc_k(M) complex
        minimize(sum_square_abs(Ck*w_wc_k))
        subject to
            real(w_wc_k'*a0) >= 1 + epsilon*norm(w_wc_k, 2);
            imag(w_wc_k'*a0) == 0;
    cvx_end

    if isempty(strfind(cvx_status, 'Solved')) %#ok<STREMP>
        warning('CVX status at angle error %.2f deg: %s', theta_error(k), cvx_status);
    else
        SINR_wc_curve(k) = output_sinr(w_wc_k, a_k);
    end

    SINR_ideal_curve(k) = output_sinr(w_ideal_k, a_k);
    SINR_mvdr_curve(k)  = output_sinr(w_mvdr_k, a_k);
end

figure;
plot(theta_error, SINR_ideal_curve, '--', 'LineWidth', 1.5);
hold on;
plot(theta_error, SINR_mvdr_curve, '-.', 'LineWidth', 1.5);
plot(theta_error, SINR_wc_curve, '-', 'LineWidth', 1.8);

xline(-Delta_theta, ':', 'Uncertainty boundary');
xline(Delta_theta, ':', 'Uncertainty boundary');

grid on;
xlabel('Steering-angle mismatch / degree');
ylabel('Output SINR / dB');
title(sprintf('Robustness comparison, uncertainty interval = \\pm%.1f^\\circ', ...
      Delta_theta));
legend('Ideal MVDR', 'Mismatched MVDR', 'WC-RAB', 'Location', 'best');
