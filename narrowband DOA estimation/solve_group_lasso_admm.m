function [X, r_history, s_history, eps_pri_history, eps_dual_history, actual_iters] = ...
    solve_group_lasso_admm(A, Y, lambda, rho, maxiter, epsilon_abs, epsilon_rel)
% Solve group LASSO using the notebook's ADMM iteration.

X = zeros(size(A, 2), size(Y, 2)); % initialize X
Z = X; % initialize Z
U = Z; % initialize dual variable
prox = @(V, alpha) max(0, 1 - alpha ./ max(vecnorm(V, 2, 2), eps)) .* V;% proximal gradient

% rho is fixed: factor the smaller system once and reuse it.
use_small_system = size(A, 2) > size(A, 1);
if use_small_system
    R = chol(A * A' + rho * eye(size(A, 1)));
else
    R = chol(A' * A + rho * eye(size(A, 2)));
    C = A' * Y;
end

% Record residuals and their stopping thresholds.
r_history = zeros(maxiter, 1);
s_history = zeros(maxiter, 1);
eps_pri_history = zeros(maxiter, 1);
eps_dual_history = zeros(maxiter, 1);
actual_iters = maxiter;

for i=1:maxiter
    Zold = Z;

    % ADMM update steps
    Q = Z - U;
    if use_small_system
        B = Y - A * Q;
        X = Q + A' * (R \ (R' \ B));
    else
        X = R \ (R' \ (C + rho * Q));
    end
    Z = prox(X + U, lambda / rho);
    U = U + X - Z;

    % Check convergence
    r_norm = norm(X - Z, 'fro');
    s_norm = norm(rho * (Z - Zold), 'fro');
    eps_pri = sqrt(numel(X)) * epsilon_abs + epsilon_rel * max(norm(X, 'fro'), norm(Z, 'fro'));
    eps_dual = sqrt(numel(X)) * epsilon_abs + epsilon_rel * norm(rho * U, 'fro');
    r_history(i) = r_norm;
    s_history(i) = s_norm;
    eps_pri_history(i) = eps_pri;
    eps_dual_history(i) = eps_dual;

    if r_norm <= eps_pri && s_norm <= eps_dual
        actual_iters = i;
        break;
    end
end

r_history = r_history(1:actual_iters);
s_history = s_history(1:actual_iters);
eps_pri_history = eps_pri_history(1:actual_iters);
eps_dual_history = eps_dual_history(1:actual_iters);
end
