/%% 基准算法4：不使用导数信息且最高使用一阶C0信息的积分 - setff加速算法
% Write by Zhou Ruihan
% Date: 2020 / 10 / 23
% Version: 0.1 Trial V.
% BASIC VERSION: (1)The int - steff algorithm without derivative information is proposed.
clear all;
clc;
%% 定义停机准则
% 键入最大迭代步数
max_step = 200;
% 键入最大函数容忍限度
max_f_step = 1e-16;
% 键入最大相邻容忍限度
max_x_step = 1e-14;
% 键入启动点
start_point = 0.041;
%% 求解过程
[root, step, choose, wList] = IntSteff(max_step, max_f_step, max_x_step, start_point);
%% 结果展示
ShowResult(root, step, choose, wList);
%% 单元函数键入
function y = f(x)
% 在此处定义函数，要求为单变元。
% 1
y = 2*x*exp(-20) - 2*exp(-20*x) + 1;
% 2
% y = 1;
% for k = 1:10
%     y = y * (x^2 + x + k);
% end
% y = (10^(-8))*(x-1)*y;
% 3
% y = exp(21000/x)/(1.11*(10^11)*x^2) - 1;
% 4
% y = log(x) + 1/x - 100;
% 5
% y = sqrt(x^4 + 8)*sin(pi/(x^2+2))^2 + x^3 / (x^4 + 1) - sqrt(6)/2 + 8/17;
end
%% 区间键入
function [a, b] = Area()
% 在此处定义根的搜索区间，两端必须符号相反。
a = -0.4;
b = 0.71;
end
%% Steffensen算法迭代流程
function [root, i, choose, wList] = IntSteff(max_step, max_f_step, max_x_step, start_point)
% 计数
i = 1;
% 取出区间端点
[a, b] = Area();
% w = start_point;
w = (a + b) / 2;
% 记录各步的root
wList = [w];
choose = 0;
% 修改区间的最大容忍度
max_x_step = max_x_step + max([abs(a), abs(b), 1])*2^(-53);
% 开始迭代
while i <= max_step
    % 函数/导数合理性检查 
    if abs(f(w)) == Inf
        fprintf('警告! 迭代无法继续进行, 因为当前解: root = %f下, 函数值或迭代数值可能趋于 +/- Inf. \n', w);
        disp('使用上一次(或初值)求出的根作为输出, 算法终止!');
        disp('请考虑使用非函数的方法进行求解.');
        root = w;
        choose = 98;
        break
    end
    % 执行迭代
    alpha = 0.5;   % 定义alpha的值
    z = f(w);
    d = w - alpha * z^2 / (f(w + z) - z);
    wnew = w + 3*(d - w)*z / (2*z - f(d));
    if abs(wnew) == Inf
        fprintf('警告！当前的迭代方向发散: %f, 使用上一次(或初值)求出的根作为输出, 算法终止! \n', w);
        disp('请检查您规定的区间范围和初始迭代点是否合理.');
        root = w;
        choose = 97;
        break
    end
    % 区间检查
    if wnew > b || wnew < a
        fprintf('警告！当前的根超出限定范围: %f, 使用上一次(或初值)求出的根作为输出, 算法终止! \n', w);
        disp('请检查您规定的区间范围是否合理.');
        root = w;
        choose = 99;
        break
    end
    % 迭代停机检查
    if (abs(wnew - w) < max_x_step)
        root = wnew;
        choose = 1;
        break
    elseif (abs(f(wnew)) < max_f_step)
        root = wnew;
        choose = 2;
        break
    end
    w = wnew;
    wList(end+1) = w;
    i = i + 1;
end
end
%% 结果展示
function ShowResult(root, step, choose, wList)
disp('您提交的 使用积分-Steff加速算法计算单元方程的根 任务已经接近完成...');
disp('算法执行完毕, 打印结果...');
if choose == 0
    fprintf('针对输入函数, 算法在第%d步完成迭代, 使用的停机准则为达到了最大的指定迭代步数. \n', step);
elseif choose == 1
    fprintf('针对输入函数, 算法在第%d步完成迭代, 使用的停机准则为相邻的迭代步长过小, 趋于收敛. \n', step);
elseif choose == 2
    fprintf('针对输入函数, 算法在第%d步完成迭代, 使用的停机准则为函数值满足最大的容忍限度要求. \n', step);
elseif choose == 99
    disp('迭代过程出现错误, 错误信息: 迭代产生的解超出了限定的区间范围.');
elseif choose == 98
    disp('迭代过程出现错误, 错误信息: 迭代过程中产生的函数值或迭代值可能趋于无穷.');
elseif choose == 97
    disp('迭代过程出现错误, 错误信息: 迭代过程中出现了迭代方向发散的情况，这使得迭代过程无法进行.')
end
fprintf('迭代求取方程的根为: %f. \n', root);
fprintf('在此根下的函数值为: %f. \n', f(root));
disp('打印 {步 - 误差} 与 {步 - log(误差)} 曲线...');
plot_x = 1:step;
plot_y = abs(wList(1:end) - root);
plot_y2 = log(plot_y);
% 曲线拟合
% 剔除-Inf项和Inf项
CHfind_Index = find(plot_y2 ~= Inf & plot_y2 ~= -Inf);
linear_fit = polyfit(plot_x(CHfind_Index), plot_y2(CHfind_Index), 1);
plot_y2_fit = polyval(linear_fit, plot_x(CHfind_Index));
k = linear_fit(1,1);
disp_k = strcat('Linear k = ', num2str(k));
disp('可能的线性收敛阶刻画：');
fprintf('x - log(error)的误差拟合直线斜率k = %f. \n', k);
% 高阶连续性拟合
plot_y3 = log(-plot_y2);
CHfind_Index_high = find(plot_y3 ~= Inf & plot_y3 ~= -Inf);
linear_fit_high = polyfit(plot_x(CHfind_Index_high), plot_y3(CHfind_Index_high), 1);
plot_y3_fit = polyval(linear_fit_high, plot_x(CHfind_Index_high));
k_high = linear_fit_high(1,1);
disp_k_high = strcat('Linear k = ', num2str(k_high));
disp('可能的非线性收敛阶刻画：');
fprintf('x - log(-log(error))的误差拟合直线斜率k = %f. \n', k_high);
% 图例显示
num_plot_x = floor(length(plot_x)/2) + 1;
num_plot_y = floor(length(plot_x)/4) + 1;
% % % % % 
figure
subplot(2, 2, 1)
plot(plot_x, plot_y, '-*', 'LineWidth', 2, 'MarkerSize', 5);
set(gca,'FontSize', 15, 'Fontname', 'Times New Roman');
xlabel('Iteration Step', 'Fontname', 'Times New Roman','FontSize', 15);
ylabel('Error of Solution', 'Fontname', 'Times New Roman','FontSize', 15);
title('Iteration Result: Step - Error', 'Fontname', 'Times New Roman','FontSize', 18);
subplot(2, 2, 2)
plot(plot_x, plot_y2, '-*', 'Color', [1, 0.5, 0], 'LineWidth', 2, 'MarkerSize', 5);
hold on
plot(plot_x(CHfind_Index), plot_y2_fit(CHfind_Index), '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MarkerSize', 3);
text(plot_x(num_plot_x), plot_y2(num_plot_y), disp_k, 'Color', [0, 0.5, 0.5], 'Fontname', 'Times New Roman','FontSize', 18);
set(gca,'FontSize', 15, 'Fontname', 'Times New Roman');
xlabel('Iteration Step', 'Fontname', 'Times New Roman','FontSize', 15);
ylabel('log(Error of Solution)', 'Fontname', 'Times New Roman','FontSize', 15);
title('Iteration Result: Step - log(Error)', 'Fontname', 'Times New Roman','FontSize', 18);
subplot(2, 2, 3)
plot(plot_x, plot_y3, '-*', 'Color', [1, 0.5, 0], 'LineWidth', 2, 'MarkerSize', 5);
hold on
plot(plot_x(CHfind_Index_high), plot_y3_fit(CHfind_Index_high), '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'MarkerSize', 3);
text(plot_x(num_plot_x), plot_y3(num_plot_y), disp_k_high, 'Color', [0, 0.5, 0.5], 'Fontname', 'Times New Roman','FontSize', 18);
set(gca,'FontSize', 15, 'Fontname', 'Times New Roman');
xlabel('Iteration Step', 'Fontname', 'Times New Roman','FontSize', 15);
ylabel('log(-log(Error of Solution))', 'Fontname', 'Times New Roman','FontSize', 15);
title('Iteration Result: Step - log(-log(Error))', 'Fontname', 'Times New Roman','FontSize', 18);
end