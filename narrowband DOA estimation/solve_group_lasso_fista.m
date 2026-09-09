function [X, res_norm, restart_iters, actual_iters] = solve_group_lasso_fista( ...
    A, Y, C, lambda, step_size, maxiter, epsilon, userestart)
% Solve group LASSO using the notebook's FISTA iteration.

X = zeros(size(A, 2), size(Y, 2));
Z = X;
t = 1;
prox = @(V, alpha) max(0, 1 - alpha ./ max(vecnorm(V, 2, 2), eps)) .* V;

% Record proximal gradient residual and restart iteration indices
res_norm = zeros(maxiter, 1);
restart_iters = [];
actual_iters = maxiter;

for i = 1:maxiter
    Zold = Z;
    Xold = X;
    gradient = A' * (A * Z) - C;
    V = Z - step_size * gradient;
    X_new = prox(V, step_size * lambda);
    t_new = (1 + sqrt(1 + 4 * t^2)) / 2;
    beta = (t - 1) / t_new;
    Z = X_new + beta * (X_new - Xold);
    X = X_new;
    t = t_new;
    grad_new = A' * (A * X_new - Y);
    G_alpha = (X_new - prox( ...
        X_new - step_size * grad_new, step_size * lambda)) / step_size;

    res_norm(i) = norm(G_alpha, 'fro');

    if res_norm(i) < epsilon
        actual_iters = i;
        res_norm = res_norm(1:i);
        break;
    end
    if userestart == true
        % m = real(trace((Zold - X_new)' * (X_new - Xold)));
        m = real(sum(conj(Zold - X_new) .* (X_new - Xold), 'all'));
        if m > 0
            t = 1;
            Z = X_new;
            restart_iters(end+1) = i; %#ok<AGROW>
        end
    end
end
if actual_iters == maxiter
    res_norm = res_norm(1:maxiter);
end
end
