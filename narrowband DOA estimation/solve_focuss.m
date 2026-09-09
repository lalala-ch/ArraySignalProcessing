function [X, fit_history, rel_history, actual_iters] = solve_focuss( ...
    A, Y, X, Gram, C, lambda, p, epsilon, maxiter, stop_threshold)
% Solve FOCUSS from the supplied initial X using the notebook's iteration.

fit_history = zeros(maxiter, 1);
rel_history = zeros(maxiter, 1);
actual_iters = maxiter;

for i=1:maxiter
    Xold = X;
    W = diag((vecnorm(Xold, 2, 2).^2 + epsilon).^(p/2 - 1));
    X = (Gram + lambda*p * W) \ C;
    fit_history(i) = norm(Y - A * X, 'fro');
    rel_history(i) = (norm(X - Xold, 'fro'))/(norm(Xold, 'fro')+eps);
    if rel_history(i) < stop_threshold
        actual_iters = i;
        break;
    end
end

fit_history = fit_history(1:actual_iters);
rel_history = rel_history(1:actual_iters);
end
