%% 基准算法1：二分法
% 
clear all;
clc;
%% 定义停机准则
% 键入最大迭代步数
max_step = 200;
% 键入最大函数容忍限度
max_f_step = 1e-16;
% 键入最大区间容忍限度
max_x_step = 1e-14;
%% 求解过程
[root, step, choose, wList] = Dichotomy(max_step, max_f_step, max_x_step);
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
% 在此处定义根的搜索区间。若无法确定根的区间，请键入相同的任意两个数字。
a = -0.4;
b = 0.71;
end
%% 二分法迭代流程
function [root, i, choose, wList] = Dichotomy(max_step, max_f_step, max_x_step)
% 计数
i = 1;
% 记录各步的root
wList = [];
choose = 0;
% 取出区间端点
[a, b] = Area();
% 修改区间的最大容忍度
max_x_step = max_x_step + max([abs(a), abs(b), 1])*2^(-53);
% 迭代过程
while i <= max_step
    w = (a + b)/2;
    wList(end+1) = w;
    % 迭代停机检查
    if (abs(b-w) < max_x_step)
        root = w;
        choose = 1;
        break
    elseif (abs(f(w)) < max_f_step)
        root = w;
        choose = 2;
        break
    end
    if f(b)*f(w) <= 0
        a = w;
    else
        b = w;
    end
    i = i + 1;
end
end
%% 结果展示
function ShowResult(root, step, choose, wList)
disp('您提交的 使用二分法计算单元方程的根 任务已经接近完成...');
disp('算法执行完毕, 打印结果...');
if choose == 0
    fprintf('针对输入函数, 算法在第%d步完成迭代, 使用的停机准则为达到了最大的指定迭代步数. \n', step);
elseif choose == 1
    fprintf('针对输入函数, 算法在第%d步完成迭代, 使用的停机准则为区间过小. \n', step);
else
    fprintf('针对输入函数, 算法在第%d步完成迭代, 使用的停机准则为函数值满足最大的容忍限度要求. \n', step);
end
fprintf('迭代求取方程的根为: %f. \n', root);
fprintf('在此根下的函数值为: %f. \n', f(root));
disp('打印 {步 - 误差} 与 {步 - log(误差)} 曲线...');
plot_x = 1:step;
plot_y = abs(wList - root);
plot_x2 = log(plot_x);
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





