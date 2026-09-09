clc;
clear;
close all;

%% =========================================================
%  1. 阵列参数
% ==========================================================

lambda = 3e-2;               % 波长
M = 12;                      % 阵元数

% 半波长均匀线阵，列向量 M×1
position = (0:M-1).' * lambda/2;

% LS 误差权重
% 注意：这是 LS 中的 alpha，不是阵元间距系数
alpha_ls = 0.5;

% 导向矢量函数
% d(theta) = exp(-j*2*pi/lambda*p*sin(theta))
steering_vector = @(theta) ...
    exp(-1j*2*pi/lambda * position * sind(theta));

%% =========================================================
%  2. 设计一：直接设计 broadside，即 theta0 = 0°
% ==========================================================

theta0_broadside = 0;

% 主瓣集合：只包含 0°
theta_m_0 = theta0_broadside;

% 旁瓣集合：
% [-90°, -20°] 取 50 点
% [ 20°,  90°] 取 50 点
theta_s_0_left  = linspace(-90, -20, 50);
theta_s_0_right = linspace( 20,  90, 50);
theta_s_0 = [theta_s_0_left, theta_s_0_right];

% 主瓣方向导向矢量
d_m_0 = steering_vector(theta_m_0);       % M×1

% 旁瓣方向导向矩阵
D_s_0 = steering_vector(theta_s_0);       % M×100

% 按书中式 (2.21) 构造 G_D 和 g_D
G_D_0 = (1-alpha_ls) * (d_m_0*d_m_0') ...
      + alpha_ls     * (D_s_0*D_s_0');

g_D_0 = (1-alpha_ls) * d_m_0;

% broadside LS 权值
w_LS_0 = G_D_0 \ g_D_0;

%% =========================================================
%  3. 设计二：直接设计一个指向 20° 的 LS beamformer
% ==========================================================

theta0_direct = 20;

% 主瓣集合：只包含 20°
theta_m_20 = theta0_direct;

% 按书中设置：
% 左侧旁瓣区域 [-90°, 0°] 取 50 点
% 右侧旁瓣区域 [ 40°, 90°] 取 50 点
theta_s_20_left  = linspace(-90, 0, 50);
theta_s_20_right = linspace( 40, 90, 50);
theta_s_20 = [theta_s_20_left, theta_s_20_right];

% 主瓣方向导向矢量
d_m_20 = steering_vector(theta_m_20);     % M×1

% 旁瓣方向导向矩阵
D_s_20 = steering_vector(theta_s_20);     % M×100

% 直接构造 20° 设计对应的 G_D 和 g_D
G_D_20 = (1-alpha_ls) * (d_m_20*d_m_20') ...
       + alpha_ls     * (D_s_20*D_s_20');

g_D_20 = (1-alpha_ls) * d_m_20;

% 直接设计得到的 20° LS 权值
w_LS_20_direct = G_D_20 \ g_D_20;

%% =========================================================
%  4. 设计三：把 broadside LS beamformer 转向 20°
% ==========================================================

theta0_steer = 20;

% 20° 方向的导向矢量
d_steer_20 = steering_vector(theta0_steer);

% 按当前符号定义：
% P(theta) = w^H d(theta)
% d(theta) 使用负指数
%
% 因此转向权值为
% w_steer = w_broadside .* d(theta0)
w_LS_20_steered = w_LS_0 .* d_steer_20;

%% =========================================================
%  5. 计算三个方向图
% ==========================================================

theta_plot = -90:0.1:90;

D_plot = steering_vector(theta_plot);     % M×Ntheta

% 三个波束响应
P_LS_0          = w_LS_0'          * D_plot;
P_LS_20_direct  = w_LS_20_direct'  * D_plot;
P_LS_20_steered = w_LS_20_steered' * D_plot;

% 分别按照各自峰值归一化
BP_LS_0 = 20*log10( ...
    abs(P_LS_0)/max(abs(P_LS_0)) + eps);

BP_LS_20_direct = 20*log10( ...
    abs(P_LS_20_direct)/max(abs(P_LS_20_direct)) + eps);

BP_LS_20_steered = 20*log10( ...
    abs(P_LS_20_steered)/max(abs(P_LS_20_steered)) + eps);

%% =========================================================
%  6. 检查三个方向图的峰值位置
% ==========================================================

[~, index_peak_0] = max(abs(P_LS_0));
[~, index_peak_direct] = max(abs(P_LS_20_direct));
[~, index_peak_steered] = max(abs(P_LS_20_steered));

theta_peak_0 = theta_plot(index_peak_0);
theta_peak_direct = theta_plot(index_peak_direct);
theta_peak_steered = theta_plot(index_peak_steered);

fprintf('Broadside LS 峰值方向：       %.2f°\n', theta_peak_0);
fprintf('直接设计 20° LS 峰值方向：   %.2f°\n', theta_peak_direct);
fprintf('Broadside 转向 20° 峰值方向：%.2f°\n', theta_peak_steered);

%% =========================================================
%  7. 绘制三个独立方向图
% ==========================================================

figure;
plot(theta_plot, BP_LS_0, 'LineWidth', 1.3);
grid on;
xlim([-90 90]);
ylim([-80 0]);

xlabel('\theta (degree)');
ylabel('Magnitude response (dB)');
title('Broadside LS-Based Fixed Beamformer');

xline(0, '--', '\theta_0=0^\circ');

figure;
plot(theta_plot, BP_LS_20_direct, 'LineWidth', 1.3);
grid on;
xlim([-90 90]);
ylim([-80 0]);

xlabel('\theta (degree)');
ylabel('Magnitude response (dB)');
title('Directly Designed 20^\circ LS-Based Fixed Beamformer');

xline(20, '--', '\theta_0=20^\circ');

figure;
plot(theta_plot, BP_LS_20_steered, 'LineWidth', 1.3);
grid on;
xlim([-90 90]);
ylim([-80 0]);

xlabel('\theta (degree)');
ylabel('Magnitude response (dB)');
title('Broadside LS Beamformer Steered to 20^\circ');

xline(20, '--', '\theta_0=20^\circ');

%% =========================================================
%  8. 对比：直接设计 20° 与 broadside 转向到 20°
% ==========================================================

figure;

plot(theta_plot, BP_LS_20_direct, ...
    'LineWidth', 1.4);
hold on;

plot(theta_plot, BP_LS_20_steered, '--', ...
    'LineWidth', 1.4);

grid on;
xlim([-90 90]);
ylim([-80 0]);

xlabel('\theta (degree)');
ylabel('Magnitude response (dB)');

title('Direct LS Design vs. Beam Steering to 20^\circ');

legend('Direct LS design at 20^\circ', ...
       'Broadside LS steered to 20^\circ', ...
       'Location', 'best');

xline(20, ':', '\theta_0=20^\circ');

%% =========================================================
%  9. 输出权值
% ==========================================================

disp(' ');
disp('Broadside LS 权值 w_LS_0：');
disp(w_LS_0);

disp('直接设计的 20° LS 权值 w_LS_20_direct：');
disp(w_LS_20_direct);

disp('由 broadside 转向得到的 20° 权值 w_LS_20_steered：');
disp(w_LS_20_steered);