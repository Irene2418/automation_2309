% 模拟出一个系统，希望逼近参数真值
% 初始化
clear; clc;
a = [1 -1.9799 0.9799 ]'; b = 9.5694*1e-4;
d = 2; na = 2; nb = 0; N = 1000; % 增加数据点数N
uk = zeros(d + nb, 1);
yk = zeros(na, 1);
kxi = sqrt(1) * randn(N, 1);  
theta = [a(2:na + 1); b];
% 计算
x1 = 1; x2 = 1; x3 = 1; x4 = 0; S = 1;
for k = 1:N
    M(k) = xor(x3, x4);
    IM(k) = xor(M(k), S);
    if IM(k) == 0
        u(k) = -1;
    else
        u(k) = 1;
    end
    S = not(S);
    x4 = x3; x3 = x2; x2 = x1; x1 = M(k);
end
% 更新
for k = 1:N % 不断更新计算Y和phi
    Phi(k, :) = [-yk; uk(d : d + nb)]';
    y(k) = Phi(k, :) * theta + kxi(k);
    for i = d + nb :-1:2
        uk(i) = uk(i - 1);
    end
    u(1) = u(k);
    for i = na :-1:2
        yk(i) = yk(i - 1);
    end
    yk(1) = y(k);
end
thetae = pinv(Phi' * Phi) * Phi' * y';

% 计算预测值
y_hat = Phi * thetae;

% 绘制实际值和预测值的曲线
figure;
plot(1:N, y, 'b', 'DisplayName', 'Actual y');
hold on;
plot(1:N, y_hat, 'r--', 'DisplayName', 'Predicted y\_hat');
xlabel('Time Step');
ylabel('Output');
legend;
title('Actual Output vs Predicted Output');
hold off;
