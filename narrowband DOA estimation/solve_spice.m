function [q, R, objective_history, rel_history, actual_iters, converged] = ...
    solve_spice(A, Y, maxiter, eta)
% Solve SPICE using the notebook's unified covariance-fitting iteration.
% A: M-by-K signal dictionary; Y: M-by-L observed snapshots.
% maxiter: positive integer iteration limit; eta: relative-change threshold.
% q: [K signal powers; M sensor noise powers]; R: fitted covariance.
% objective_history includes the initial value and each updated iterate.
% rel_history contains one relative power change per iteration.
% converged is true only when the relative change falls below eta.

[M, L] = size(Y);
B = [A eye(M)];
num_powers = size(B, 2);
R_SCM = Y*Y'/L;
R_SCM = (R_SCM + R_SCM')/2;

% Positive periodogram initialization, one power per dictionary column.
b_norm2 = sum(abs(B).^2, 1).';
power_scale = real(trace(R_SCM))/M;
assert(power_scale > 0, 'SPICE requires nonzero observed data.');
power_floor = 1e-8*power_scale;
q = real(sum(conj(B).*(R_SCM*B), 1)).' ./ (b_norm2.^2);
q = max(q, power_floor);

is_singular = rank(R_SCM) < M;
if is_singular
    V = R_SCM;
    W = b_norm2;
else
    V = sqrtm(R_SCM);
    V = (V + V')/2;
    W = real(sum(conj(B).*(R_SCM\B), 1)).';
end
sqrt_W = sqrt(W); % fixed throughout the iteration

objective_history = nan(maxiter + 1, 1);
rel_history = nan(maxiter, 1);
converged = false;
for i = 1:maxiter
    qold = q;
    Q = spdiags(q, 0, num_powers, num_powers);
    R = B*Q*B';
    R = (R + R')/2;
    U = R\V;
    C = Q*(B'*U);
    objective_history(i) = real(trace(V'*U) + W'*qold);

    q = vecnorm(C, 2, 2)./sqrt_W;
    rel_history(i) = norm(q-qold, 2)/max(norm(qold, 2), eps);
    if rel_history(i) < eta
        converged = true;
        break;
    end
end
actual_iters = i;
rel_history = rel_history(1:actual_iters);

% Recompute the model and objective for the final, updated power vector.
Q = spdiags(q, 0, num_powers, num_powers);
R = B*Q*B';
R = (R + R')/2;
U = R\V;
objective_history(actual_iters + 1) = real(trace(V'*U) + W'*q);
objective_history = objective_history(1:actual_iters + 1);
end
