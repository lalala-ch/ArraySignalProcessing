function [gamma, sigma2, mu, log_evidence_history, rel_history, ...
    actual_iters, converged] = solve_sbl_em(A, Y, maxiter, eta)
% Solve multisnapshot EM-SBL with a common angular prior and noise variance.
% A: M-by-N steering dictionary; Y: M-by-L observed snapshots.
% maxiter: positive integer iteration limit; eta: relative-change threshold.
% gamma: N angular prior variances; sigma2: common sensor noise variance.
% mu: N-by-L posterior means computed from the final hyperparameters.
% log_evidence_history includes the initial and each updated iterate.
% rel_history records the larger relative change in gamma and sigma2.

[M, L] = size(Y);
num_grid = size(A, 2);
assert(M > 0 && L > 0 && num_grid > 0, ...
    'A and Y must contain sensors, grid points, and snapshots.');
assert(size(A, 1) == M, 'A and Y must have the same number of sensor rows.');
assert(isscalar(maxiter) && maxiter >= 1 && maxiter == floor(maxiter), ...
    'maxiter must be a positive integer.');
assert(isscalar(eta) && eta > 0, 'eta must be positive.');

% Initialize the model covariance at the scale of the observed power.
P_mean = real(sum(abs(Y(:)).^2)/(M*L));
assert(P_mean > 0, 'EM-SBL requires nonzero observed data.');
power_floor = 1e-12*max(P_mean, 1e-12);
ka = real(sum(abs(A(:)).^2)/M);
assert(ka > 0, 'The steering dictionary A must be nonzero.');
sigma2 = max(P_mean/10, power_floor);
gamma_0 = max(0.9*P_mean/ka, power_floor);
gamma = gamma_0*ones(num_grid, 1);

log_evidence_history = nan(maxiter + 1, 1);
rel_history = nan(maxiter, 1);
converged = false;
for iter = 1:maxiter
    % E-step: all snapshots share Gamma and the posterior covariance.
    % Compute only diag(Sigma), avoiding a num_grid-by-num_grid matrix.
    C = (A .* gamma.')*A' + sigma2*eye(M);
    C = (C + C')/2;
    CY = C\Y;
    CA = C\A;
    mu = gamma .* (A'*CY);
    aCa = real(sum(conj(A).*CA, 1)).';
    sigma_diag = max(real(gamma - gamma.^2 .* aCa), 0);

    R = chol(C);
    logdet_C = 2*sum(log(diag(R)));
    log_evidence_history(iter) = -M*L*log(pi) - L*logdet_C ...
        - real(sum(conj(Y).*CY, 'all'));

    % M-step: both updates use the same old posterior.
    gamma_new = max(sum(abs(mu).^2, 2)/L + sigma_diag, power_floor);
    % tr(A*Sigma*A') = sigma2*sum(gamma.*aCa).
    trace_ASigmaAH = sigma2*sum(gamma .* aCa);
    sigma2_new = max(real((norm(Y-A*mu, 'fro')^2 + ...
        L*trace_ASigmaAH)/(M*L)), power_floor);

    rel_gamma = norm(gamma_new-gamma, 2)/max(norm(gamma, 2), power_floor);
    rel_sigma = abs(sigma2_new-sigma2)/max(sigma2, power_floor);
    rel_history(iter) = max(rel_gamma, rel_sigma);
    gamma = gamma_new;
    sigma2 = sigma2_new;
    if rel_history(iter) < eta
        converged = true;
        break;
    end
end

actual_iters = iter;
rel_history = rel_history(1:actual_iters);

% Return a posterior and evidence value consistent with the final parameters.
C = (A .* gamma.')*A' + sigma2*eye(M);
C = (C + C')/2;
CY = C\Y;
mu = gamma .* (A'*CY);
R = chol(C);
logdet_C = 2*sum(log(diag(R)));
log_evidence_history(actual_iters + 1) = -M*L*log(pi) ...
    - L*logdet_C - real(sum(conj(Y).*CY, 'all'));
log_evidence_history = log_evidence_history(1:actual_iters + 1);
end
