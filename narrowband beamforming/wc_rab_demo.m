clc;
clear;
close all;
rng(7);

%% Worst-case robust adaptive beamforming (WC-RAB)
% 实现问题：
%   min_w  w^H Rhat w
%   s.t.   Re{w^H a0} >= 1 + epsilon*||w||_2
%          Im{w^H a0} = 0
%
% 本代码不依赖 CVX。
% 利用等效对角加载形式：
%   q(mu) = (Rhat + mu I)^(-1) a0
%   mu * ||q(mu)||_2 = epsilon
%   w = q(mu)/(a0^H q(mu) - epsilon||q(mu)||_2)

%% 1. 阵列与信号参数
M = 10;                         % 阵元数
N = 300;                        % 快拍数
m = (0:M-1).';                  % 阵元编号
d_over_lambda = 0.5;            % 阵元间距/波长

theta0 = 0;                     % 假定的目标方向
theta_s_true = 2;               % 真实目标方向
theta_i = [-30, 40];            % 两个干扰方向

SNR_dB = 10;                    % 单阵元输入信噪比
INR_dB = [30, 25];              % 两个干扰的单阵元干噪比
sigma2 = 1;                     % 单阵元噪声功率

Ps = sigma2 * 10^(SNR_dB/10);
Pi = sigma2 * 10.^(INR_dB/10);

% ULA 导向矢量；角度相对于阵列 broadside
steer = @(theta) exp(-1j*2*pi*d_over_lambda*m*sind(theta));

a0 = steer(theta0);
a_true = steer(theta_s_true);
a_i1 = steer(theta_i(1));
a_i2 = steer(theta_i(2));

%% 2. 根据角度不确定区间设置 epsilon
% 假设真实目标方向位于 [theta0-Delta_theta, theta0+Delta_theta]
Delta_theta = 2.5;              % 角度不确定半宽，单位：度
theta_unc = linspace(theta0-Delta_theta, theta0+Delta_theta, 2001);
A_unc = steer(theta_unc);

% 用一个以 a0 为中心的二范数球包住整个角度区间
epsilon = max(vecnorm(A_unc - a0, 2, 1));

fprintf('||a0||_2 = %.6f\n', norm(a0));
fprintf('epsilon  = %.6f\n', epsilon);
fprintf('实际导向矢量误差 ||a_true-a0||_2 = %.6f\n', norm(a_true-a0));

if epsilon >= norm(a0)
    error(['epsilon 必须小于 ||a0||，否则鲁棒约束不可行。', ...
           '请减小 Delta_theta，或重新检查导向矢量的归一化方式。']);
end

if norm(a_true-a0) <= epsilon
    fprintf('真实导向矢量位于所设置的不确定球内。\n\n');
else
    warning('真实导向矢量不在所设置的不确定球内。');
end

%% 3. 生成阵列接收数据
% 目标信号：QPSK，单位平均功率
data = randi([0, 3], 1, N);
s_unit = exp(1j*(pi/4 + pi/2*data));
s = sqrt(Ps) * s_unit;

% 两路干扰：独立复高斯信号
u1 = (randn(1,N) + 1j*randn(1,N))/sqrt(2);
u2 = (randn(1,N) + 1j*randn(1,N))/sqrt(2);
i1 = sqrt(Pi(1)) * u1;
i2 = sqrt(Pi(2)) * u2;

% 阵元噪声：空间白、时间白复高斯噪声
noise_unit = (randn(M,N) + 1j*randn(M,N))/sqrt(2);
noise = sqrt(sigma2) * noise_unit;

% M x N 接收数据矩阵
X = a_true*s + a_i1*i1 + a_i2*i2 + noise;

% 样本协方差矩阵
Rhat = X*X'/N;
Rhat = (Rhat + Rhat')/2;

% 极小的数值正则项，避免有限精度下病态
reg = 1e-8 * real(trace(Rhat))/M;
Rhat = Rhat + reg*eye(M);

%% 4. 三种波束形成器

% 4.1 理想 MVDR：设计时知道真实目标方向，仅作为性能参考
q = Rhat \ a_true;
w_ideal = q/(a_true'*q);

% 4.2 失配 MVDR：错误地使用 theta0
q = Rhat \ a0;
w_mvdr = q/(a0'*q);

% 4.3 Worst-case 鲁棒自适应波束形成器
[w_wc, mu_wc] = wc_rab_bisection(Rhat, a0, epsilon);

%% 5. 检查约束
nominal_response = w_wc' * a0;
worst_case_lower_bound = real(nominal_response) - epsilon*norm(w_wc);

fprintf('WC-RAB 等效对角加载参数 mu = %.6e\n', mu_wc);
fprintf('Re{w_wc^H a0} = %.9f\n', real(nominal_response));
fprintf('Im{w_wc^H a0} = %.3e\n', imag(nominal_response));
fprintf('Re{w_wc^H a0} - epsilon||w_wc|| = %.9f\n\n', ...
        worst_case_lower_bound);

%% 6. 理论输出 SINR
Rin = Pi(1)*(a_i1*a_i1') + ...
      Pi(2)*(a_i2*a_i2') + ...
      sigma2*eye(M);

output_sinr = @(w,at) ...
    10*log10(Ps*abs(w'*at)^2 / real(w'*Rin*w));

SINR_ideal = output_sinr(w_ideal, a_true);
SINR_mvdr  = output_sinr(w_mvdr,  a_true);
SINR_wc    = output_sinr(w_wc,    a_true);

fprintf('理想 MVDR 输出 SINR   = %8.3f dB\n', SINR_ideal);
fprintf('失配 MVDR 输出 SINR   = %8.3f dB\n', SINR_mvdr);
fprintf('WC-RAB 输出 SINR       = %8.3f dB\n', SINR_wc);
fprintf('|w_mvdr^H a_true|      = %.6f\n', abs(w_mvdr'*a_true));
fprintf('|w_wc^H a_true|        = %.6f\n', abs(w_wc'*a_true));

%% 7. 绘制波束图
theta_scan = -90:0.1:90;
A_scan = steer(theta_scan);

B_ideal = abs(w_ideal' * A_scan);
B_mvdr  = abs(w_mvdr'  * A_scan);
B_wc    = abs(w_wc'    * A_scan);

B_ideal_dB = 20*log10(B_ideal/max(B_ideal) + 1e-12);
B_mvdr_dB  = 20*log10(B_mvdr/max(B_mvdr)   + 1e-12);
B_wc_dB    = 20*log10(B_wc/max(B_wc)       + 1e-12);

figure;
plot(theta_scan, B_ideal_dB, '--', 'LineWidth', 1.4);
hold on;
plot(theta_scan, B_mvdr_dB,  '-.', 'LineWidth', 1.4);
plot(theta_scan, B_wc_dB,    '-',  'LineWidth', 1.8);

xline(theta0, ':', 'Presumed SOI', 'LineWidth', 1.1);
xline(theta_s_true, ':', 'True SOI', 'LineWidth', 1.1);
xline(theta_i(1), ':', 'Interference 1', 'LineWidth', 1.1);
xline(theta_i(2), ':', 'Interference 2', 'LineWidth', 1.1);

grid on;
xlim([-90, 90]);
ylim([-80, 5]);
xlabel('Direction / degree');
ylabel('Normalized response / dB');
title(sprintf('Beam patterns: true SOI = %.1f^\\circ, presumed SOI = %.1f^\\circ', ...
      theta_s_true, theta0));
legend('Ideal MVDR', 'Mismatched MVDR', 'WC-RAB', ...
       'Location', 'southoutside', 'NumColumns', 3);

%% 8. 扫描目标角度失配，比较输出 SINR
theta_error = -5:0.25:5;
SINR_ideal_curve = zeros(size(theta_error));
SINR_mvdr_curve  = zeros(size(theta_error));
SINR_wc_curve    = zeros(size(theta_error));

for k = 1:numel(theta_error)
    theta_k = theta0 + theta_error(k);
    a_k = steer(theta_k);

    % 使用同一组源信号、干扰和噪声，只改变目标真实方向
    Xk = a_k*s + a_i1*i1 + a_i2*i2 + noise;

    Rk = Xk*Xk'/N;
    Rk = (Rk + Rk')/2;
    Rk = Rk + 1e-8*real(trace(Rk))/M*eye(M);

    % 理想 MVDR
    q = Rk \ a_k;
    w_ideal_k = q/(a_k'*q);

    % 失配 MVDR
    q = Rk \ a0;
    w_mvdr_k = q/(a0'*q);

    % WC-RAB
    w_wc_k = wc_rab_bisection(Rk, a0, epsilon);

    SINR_ideal_curve(k) = output_sinr(w_ideal_k, a_k);
    SINR_mvdr_curve(k)  = output_sinr(w_mvdr_k,  a_k);
    SINR_wc_curve(k)    = output_sinr(w_wc_k,    a_k);
end

figure;
plot(theta_error, SINR_ideal_curve, '--', 'LineWidth', 1.5);
hold on;
plot(theta_error, SINR_mvdr_curve,  '-.', 'LineWidth', 1.5);
plot(theta_error, SINR_wc_curve,    '-',  'LineWidth', 1.8);

xline(-Delta_theta, ':', 'Uncertainty boundary', 'LineWidth', 1.1);
xline( Delta_theta, ':', 'Uncertainty boundary', 'LineWidth', 1.1);

grid on;
xlabel('Steering-angle mismatch / degree');
ylabel('Output SINR / dB');
title(sprintf('Robustness versus steering mismatch, uncertainty interval = \\pm%.1f^\\circ', ...
      Delta_theta));
legend('Ideal MVDR', 'Mismatched MVDR', 'WC-RAB', ...
       'Location', 'best');

%% 9. 可选：用 CVX 直接求解并与二分法结果核对
% 安装 CVX 后取消下面代码的注释。
%
% cvx_begin quiet
%     variable w_cvx(M) complex
%     minimize( real(quad_form(w_cvx, Rhat)) )
%     subject to
%         real(w_cvx' * a0) >= 1 + epsilon*norm(w_cvx, 2);
%         imag(w_cvx' * a0) == 0;
% cvx_end
%
% phase_align = exp(-1j*angle(w_cvx' * w_wc));
% relative_error = norm(w_cvx - phase_align*w_wc)/norm(w_cvx);
% fprintf('CVX 与二分法解的相对误差 = %.3e\n', relative_error);


%% 局部函数：不用 CVX 求解 WC-RAB
function [w, mu] = wc_rab_bisection(R, a0, epsilon)
%WC_RAB_BISECTION 求解球形导向矢量不确定集下的 WC-RAB
%
% 输入：
%   R       : Hermitian 正定样本协方差矩阵
%   a0      : 假定导向矢量
%   epsilon : 不确定球半径，必须满足 0 <= epsilon < ||a0||
%
% 输出：
%   w       : 鲁棒波束形成权向量
%   mu      : 等效对角加载参数
%
% 原理：
%   q(mu) = (R + mu I)^(-1)a0
%   求解 mu||q(mu)|| = epsilon
%   然后利用有效约束归一化：
%   w = q/(a0^Hq - epsilon||q||)

    R = (R + R')/2;
    M = length(a0);
    I = eye(M);

    if epsilon < 0
        error('epsilon 不能为负数。');
    end

    if epsilon >= norm(a0)
        error('必须满足 epsilon < ||a0||，否则鲁棒约束不可行。');
    end

    % epsilon = 0 时退化为普通 MVDR
    if epsilon == 0
        q = R \ a0;
        w = q/(a0'*q);
        mu = 0;
        return;
    end

    f = @(x) x*norm((R + x*I)\a0, 2) - epsilon;

    % 搜索包含根的区间 [mu_low, mu_high]
    mu_low = 0;
    mu_high = max(real(trace(R))/M, 1e-8);

    count = 0;
    while f(mu_high) < 0
        mu_high = 2*mu_high;
        count = count + 1;

        if count > 200 || ~isfinite(mu_high)
            error('未能找到二分搜索上界，请检查 R、a0 和 epsilon。');
        end
    end

    % 二分搜索
    max_iter = 100;
    rel_tol = 1e-12;

    for iter = 1:max_iter
        mu_mid = (mu_low + mu_high)/2;

        if f(mu_mid) < 0
            mu_low = mu_mid;
        else
            mu_high = mu_mid;
        end

        if (mu_high - mu_low) <= ...
                rel_tol*max(1, abs(mu_low) + abs(mu_high))
            break;
        end
    end

    mu = (mu_low + mu_high)/2;
    q = (R + mu*I) \ a0;

    denominator = real(a0'*q) - epsilon*norm(q, 2);

    if denominator <= 0
        error('归一化分母非正，请检查输入参数或矩阵条件数。');
    end

    w = q/denominator;

    % 消除数值误差造成的极小整体相位
    w = w * exp(1j*angle(w'*a0));
end
