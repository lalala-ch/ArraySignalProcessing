clc;
clear;
close all;
rng(42);

%% ============================================================
% Gu & Leshem, IEEE TSP 2012
% Robust adaptive beamforming based on INCM reconstruction
% and steering-vector estimation
%
% This script implements:
%   (6)  Capon spatial spectrum from the total sample covariance
%   (7)  interference-plus-noise covariance reconstruction
%   (11) steering-vector mismatch estimation using CVX
%   (12) corrected steering vector
%   (13) final reconstruction-plus-estimation beamformer
%
% It also compares:
%   1. Capon spectrum from total sample covariance Rhat_x
%   2. Capon-type spectrum from reconstructed Rin_rec
%   3. Capon-type spectrum from exact Rin_true
%
% Requirements:
%   MATLAB + CVX (run cvx_setup before executing this script)
%% ============================================================

%% 0. Check CVX
if exist('cvx_begin','file') ~= 2
    error(['CVX was not found. Install CVX, run cvx_setup, ', ...
           'and then execute this script again.']);
end

%% 1. Array model: half-wavelength ULA
M = 12;
m = (0:M-1).';

theta_scan = -90:0.1:90;
A_scan = exp(-1j*pi*m*sind(theta_scan));

steer = @(theta) exp(-1j*pi*m*sind(theta));

%% 2. DOAs and desired-signal angular sector
theta_s_nominal = 5;      % presumed/nominal desired DOA
theta_s_true    = 8;      % actual desired DOA: 3-degree mismatch

theta_i1 = -50;
theta_i2 = -20;

a_bar  = steer(theta_s_nominal);
a_true = steer(theta_s_true);
a_i1   = steer(theta_i1);
a_i2   = steer(theta_i2);

% Prior desired-signal sector.
% It must contain the desired signal and exclude the interferers.
Theta = [0, 10];

idx_out_sector = theta_scan < Theta(1) | theta_scan > Theta(2);

%% 3. Source powers and snapshots
sigma2 = 1;

SNR_dB  = 20;
INR1_dB = 30;
INR2_dB = 30;

Ps  = sigma2 * 10^(SNR_dB/10);
Pi1 = sigma2 * 10^(INR1_dB/10);
Pi2 = sigma2 * 10^(INR2_dB/10);

N = 500;
% Try N = 30 to inspect the finite-snapshot case.

%% 4. Generate narrowband complex Gaussian data
s  = (randn(1,N) + 1j*randn(1,N))/sqrt(2);
i1 = (randn(1,N) + 1j*randn(1,N))/sqrt(2);
i2 = (randn(1,N) + 1j*randn(1,N))/sqrt(2);

noise = sqrt(sigma2/2) * ...
        (randn(M,N) + 1j*randn(M,N));

X = sqrt(Ps)  * a_true * s ...
  + sqrt(Pi1) * a_i1   * i1 ...
  + sqrt(Pi2) * a_i2   * i2 ...
  + noise;

%% 5. Exact covariance matrices (simulation ground truth only)
Rin_true = Pi1*(a_i1*a_i1') ...
         + Pi2*(a_i2*a_i2') ...
         + sigma2*eye(M);

Rx_true = Ps*(a_true*a_true') + Rin_true;

%% 6. Total sample covariance
Rhat_x = X*X'/N;
Rhat_x = (Rhat_x + Rhat_x')/2;

% Tiny loading for stable spectral inversion only.
delta_x = 1e-6 * real(trace(Rhat_x))/M;
Rhat_x_loaded = Rhat_x + delta_x*eye(M);

%% 7. Equation (6): Capon spectrum from the total sample covariance
P_total = capon_spectrum(Rhat_x_loaded, A_scan);

%% 8. Equation (7): reconstruct the interference-plus-noise covariance
%
% Rin_rec ≈ sum_{theta_q outside Theta}
%           P_total(theta_q)d(theta_q)d^H(theta_q)Delta_theta
%
dtheta_rad = deg2rad(theta_scan(2)-theta_scan(1));

Rin_rec = zeros(M,M);
out_indices = find(idx_out_sector);

for kk = 1:numel(out_indices)
    q = out_indices(kk);
    dq = A_scan(:,q);

    Rin_rec = Rin_rec ...
        + P_total(q)*(dq*dq')*dtheta_rad;
end

Rin_rec = (Rin_rec + Rin_rec')/2;

% Numerical eigenvalue flooring is needed because CVX requires
% Hermitian positive-(semi)definite quadratic-form matrices.
Rin_rec = make_hpd(Rin_rec, 1e-8);

%% 9. Equation (8): reconstructed INCM + nominal steering vector
w_rec_nominal = mvdr_weight(Rin_rec, a_bar);

%% 10. Equation (11): CVX steering-vector estimation
%
% min_e  (a_bar+e)^H Rin_rec^{-1}(a_bar+e)
%
% s.t.   a_bar^H e = 0
%
%        (a_bar+e)^H Rin_rec(a_bar+e)
%        <= a_bar^H Rin_rec a_bar
%
Rin_rec_inv = Rin_rec \ eye(M);
Rin_rec_inv = (Rin_rec_inv + Rin_rec_inv')/2;
Rin_rec_inv = make_hpd(Rin_rec_inv, 1e-10);

constraint_bound = real(a_bar' * Rin_rec * a_bar);

cvx_begin 
    % cvx_precision high
    % cvx_solver sedumi
    variable e_perp(M) complex

    minimize( quad_form(a_bar + e_perp, Rin_rec_inv) )

    subject to
        % Complex orthogonality constraint a_bar^H e_perp = 0
        real(a_bar' * e_perp) == 0;
        imag(a_bar' * e_perp) == 0;
        % a_bar' * e_perp == 0;
        % norm(a_bar + e_perp) == sqrt(M);
        % Interference-region weighted-correlation constraint
        quad_form(a_bar + e_perp, Rin_rec) ...
            <= constraint_bound;
cvx_end

if ~(contains(cvx_status,'Solved') || ...
     contains(cvx_status,'Inaccurate/Solved'))
    error('CVX failed. Status: %s', cvx_status);
end

%% 11. Equations (12) and (13): corrected steering vector and final weight
a_est = a_bar + e_perp;

w_rec_est = mvdr_weight(Rin_rec, a_est);

%% 12. Reference beamformers for comparison
% Oracle MVDR: exact Rin and exact desired steering vector
w_oracle = mvdr_weight(Rin_true, a_true);

% Conventional SMI/MPDR: total sample covariance + nominal steering
w_smi_nominal = mvdr_weight(Rhat_x_loaded, a_bar);

% Exact MPDR: exact total covariance + exact desired steering
% It should coincide with oracle MVDR.
w_mpdr_exact = mvdr_weight(Rx_true, a_true);

%% 13. CVX and steering-vector diagnostics
orthogonality_residual = abs(a_bar' * e_perp);

constraint_lhs = real(a_est' * Rin_rec * a_est);
constraint_rhs = real(a_bar' * Rin_rec * a_bar);
constraint_residual = constraint_lhs - constraint_rhs;

rho_nominal = abs(a_true' * a_bar) ...
            / (norm(a_true)*norm(a_bar));

rho_est = abs(a_true' * a_est) ...
        / (norm(a_true)*norm(a_est));

% Project the estimated vector onto the ideal ULA manifold for diagnosis.
% The paper estimates a complex vector, not an angle.
manifold_similarity = abs(sum(conj(A_scan).*a_est,1)) ...
                    ./ (vecnorm(A_scan,2,1)*norm(a_est));

[~,idx_equiv] = max(manifold_similarity);
theta_est_equiv = theta_scan(idx_equiv);

relative_mvdr_mpdr_error = ...
    norm(w_oracle-w_mpdr_exact)/norm(w_oracle);

fprintf('CVX status                         = %s\n', cvx_status);
fprintf('CVX objective value                = %.6e\n', cvx_optval);
fprintf('||e_perp||_2                       = %.6f\n', norm(e_perp));
fprintf('|a_bar^H e_perp|                   = %.3e\n', ...
    orthogonality_residual);
fprintf('Quadratic-constraint residual      = %.3e\n', ...
    constraint_residual);
fprintf('Nominal/true steering similarity   = %.6f\n', rho_nominal);
fprintf('Estimated/true steering similarity = %.6f\n', rho_est);
fprintf('Equivalent ideal-manifold angle    = %.2f degree\n', ...
    theta_est_equiv);
fprintf('Exact MVDR/MPDR weight difference  = %.3e\n\n', ...
    relative_mvdr_mpdr_error);

%% 14. Output SINR comparison
SINR_oracle      = output_sinr(w_oracle,      Ps, a_true, Rin_true);
SINR_smi_nominal = output_sinr(w_smi_nominal, Ps, a_true, Rin_true);
SINR_rec_nominal = output_sinr(w_rec_nominal, Ps, a_true, Rin_true);
SINR_rec_est     = output_sinr(w_rec_est,     Ps, a_true, Rin_true);

fprintf('Output SINR:\n');
fprintf('  Oracle MVDR                           = %8.3f dB\n', ...
    SINR_oracle);
fprintf('  SMI/MPDR + nominal steering           = %8.3f dB\n', ...
    SINR_smi_nominal);
fprintf('  Reconstructed Rin + nominal steering  = %8.3f dB\n', ...
    SINR_rec_nominal);
fprintf('  Reconstructed Rin + CVX estimate      = %8.3f dB\n', ...
    SINR_rec_est);

%% 15. Requested three spectra
% (1) Total-data Capon spectrum
% P_total was calculated in Section 7.

% (2) Capon-type spectrum from reconstructed Rin
P_rec_type = capon_spectrum(Rin_rec, A_scan);

% (3) Capon-type spectrum from exact Rin
P_true_type = capon_spectrum(Rin_true, A_scan);

%% 16. Figure 1: individually normalized spectral shapes
P_total_norm_dB = 10*log10(P_total/max(P_total) + eps);
P_rec_norm_dB   = 10*log10(P_rec_type/max(P_rec_type) + eps);
P_true_norm_dB  = 10*log10(P_true_type/max(P_true_type) + eps);

figure;
plot(theta_scan, P_total_norm_dB, 'LineWidth', 1.5);
hold on;
plot(theta_scan, P_rec_norm_dB, '--', 'LineWidth', 1.7);
plot(theta_scan, P_true_norm_dB, '-.', 'LineWidth', 1.7);

xline(theta_s_nominal, ':', 'Nominal SOI');
xline(theta_s_true,    ':', 'True SOI');
xline(theta_i1,        ':', 'Interference 1');
xline(theta_i2,        ':', 'Interference 2');
xline(Theta(1),        ':', '\Theta lower');
xline(Theta(2),        ':', '\Theta upper');

grid on;
xlim([-90,90]);
ylim([-70,5]);
xlabel('DOA / degree');
ylabel('Individually normalized spectrum / dB');
title('Capon spectrum and INCM Capon-type spectra');
legend('Capon spectrum from total sample covariance \hat{R}_x', ...
       'Capon-type spectrum from reconstructed \tilde{R}_{i+n}', ...
       'Capon-type spectrum from exact R_{i+n}', ...
       'Location','best');

%% 17. Figure 2: absolute spectral comparison
%
% The angular integral can introduce an overall scale into Rin_rec.
% A positive scale does not alter the solution of equation (11) or
% the final MVDR weight, but it scales the Capon-type spectrum.
%
% Ground-truth trace matching is therefore used only for this
% simulation display. It is not part of the practical algorithm.
trace_scale = real(trace(Rin_true))/real(trace(Rin_rec));
Rin_rec_trace_matched = trace_scale * Rin_rec;

P_rec_abs = capon_spectrum(Rin_rec_trace_matched, A_scan);

noise_floor = sigma2/M;

P_total_abs_dB = 10*log10(P_total/noise_floor + eps);
P_rec_abs_dB   = 10*log10(P_rec_abs/noise_floor + eps);
P_true_abs_dB  = 10*log10(P_true_type/noise_floor + eps);

figure;
plot(theta_scan, P_total_abs_dB, 'LineWidth', 1.5);
hold on;
plot(theta_scan, P_rec_abs_dB, '--', 'LineWidth', 1.7);
plot(theta_scan, P_true_abs_dB, '-.', 'LineWidth', 1.7);

xline(theta_s_nominal, ':', 'Nominal SOI');
xline(theta_s_true,    ':', 'True SOI');
xline(theta_i1,        ':', 'Interference 1');
xline(theta_i2,        ':', 'Interference 2');
yline(0, ':', 'White-noise Capon floor');

grid on;
xlim([-90,90]);
xlabel('DOA / degree');
ylabel('Spectrum relative to \sigma_n^2/M / dB');
title('Absolute spectral comparison');
legend('Capon spectrum from \hat{R}_x', ...
       'Capon-type spectrum from trace-matched \tilde{R}_{i+n}', ...
       'Capon-type spectrum from exact R_{i+n}', ...
       'Location','best');

%% 18. Figure 3: beam-pattern comparison
B_oracle      = absolute_beampattern(w_oracle,      A_scan);
B_smi_nominal = absolute_beampattern(w_smi_nominal, A_scan);
B_rec_nominal = absolute_beampattern(w_rec_nominal, A_scan);
B_rec_est     = absolute_beampattern(w_rec_est,     A_scan);

figure;
plot(theta_scan, B_oracle, 'LineWidth', 1.8);
hold on;
plot(theta_scan, B_smi_nominal, '--', 'LineWidth', 1.4);
plot(theta_scan, B_rec_nominal, '-.', 'LineWidth', 1.5);
plot(theta_scan, B_rec_est, ':', 'LineWidth', 2.0);

xline(theta_s_nominal, ':', 'Nominal SOI');
xline(theta_s_true,    ':', 'True SOI');
xline(theta_i1,        ':', 'Interference 1');
xline(theta_i2,        ':', 'Interference 2');
yline(0, ':', 'Unity response');

grid on;
xlim([-90,90]);
ylim([-80,10]);
xlabel('DOA / degree');
ylabel('|w^H d(\theta)| / dB');
title('Beam-pattern comparison');
legend('Oracle MVDR', ...
       'SMI/MPDR + nominal steering', ...
       'Reconstructed Rin + nominal steering', ...
       'Reconstructed Rin + CVX steering estimate', ...
       'Location','best');

%% 19. Figure 4: steering-vector similarity diagnostics
similarity_nominal = abs(sum(conj(A_scan).*a_bar,1)) ...
                   ./ (vecnorm(A_scan,2,1)*norm(a_bar));

similarity_true = abs(sum(conj(A_scan).*a_true,1)) ...
                ./ (vecnorm(A_scan,2,1)*norm(a_true));

figure;
plot(theta_scan, similarity_nominal, 'LineWidth', 1.4);
hold on;
plot(theta_scan, similarity_true, '--', 'LineWidth', 1.4);
plot(theta_scan, manifold_similarity, '-.', 'LineWidth', 1.7);

xline(theta_s_nominal, ':', 'Nominal SOI');
xline(theta_s_true,    ':', 'True SOI');
xline(theta_est_equiv, ':', 'Equivalent estimated angle');

grid on;
xlim([Theta(1)-5, Theta(2)+5]);
ylim([0,1.05]);
xlabel('DOA / degree');
ylabel('Normalized steering-vector similarity');
title('Projection of steering vectors onto the ideal ULA manifold');
legend('Nominal steering vector', ...
       'True steering vector', ...
       'CVX-corrected steering vector', ...
       'Location','best');

%% ============================================================
% Local functions
%% ============================================================
function P = capon_spectrum(R, A)
    R = (R + R')/2;
    Z = R \ A;

    denominator = real(sum(conj(A).*Z,1));
    denominator = max(denominator, eps);

    P = 1 ./ denominator;
end

function w = mvdr_weight(R, a)
    R = (R + R')/2;
    q = R \ a;

    denominator = a' * q;
    if abs(denominator) < eps
        error('The MVDR normalization denominator is too small.');
    end

    w = q / denominator;
end

function sinr_dB = output_sinr(w, Ps, a_true, Rin_true)
    signal_power = Ps * abs(w' * a_true)^2;
    interference_noise_power = real(w' * Rin_true * w);

    sinr_dB = 10*log10(signal_power/interference_noise_power);
end

function B_dB = absolute_beampattern(w, A)
    B = abs(w' * A);
    B_dB = 20*log10(B + 1e-12);
end

function R_hpd = make_hpd(R, relative_floor)
%MAKE_HPD Force a numerically Hermitian matrix to be positive definite.

    R = (R + R')/2;

    [V,D] = eig(R,'vector');
    D = real(D);

    max_eig = max(D);
    if max_eig <= 0
        error('The matrix has no positive eigenvalue.');
    end

    eig_floor = relative_floor * max_eig;
    D = max(D,eig_floor);

    R_hpd = V*diag(D)*V';
    R_hpd = (R_hpd + R_hpd')/2;
end
