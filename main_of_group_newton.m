%% 本文实现方程组形式的牛顿迭代法。
% 本程序实现基本牛顿迭代法.
% Write by Zhou Ruihan
clear all;
clc;
%% 定义停机准则
% 键入最大迭代步数
max_step = 1000;
% 键入最大函数容忍限度
max_f_step = 1e-16;
% 键入最大相邻容忍限度
max_x_step = 1e-16;
%% 求解过程
[root, step, choose, wList] = Newton(max_step, max_f_step, max_x_step);
%% 结果展示
ShowResult(root, step, choose, wList);

%% 单元函数键入
function [y1, y2, yd] = f()
% 在此处定义函数，要求为多变元。
syms x1 x2 x3
y1 = exp(x1^2 + x2^2) - exp(2);
y2 = exp(x1^2 - x2^2) - 1;
% 计算Jacobi矩阵
yd = jacobian([y1; y2], [x1 x2]);
end
%% 计算函数值与Jacobi矩阵的值
function [yn, ydn] = f_calc(x_numeric, y1, y2, yd)
% 在此处计算函数值
syms x1 x2 x3
y1n = double(subs(y1, [x1 x2], x_numeric));
y2n = double(subs(y2, [x1 x2], x_numeric));
yn = [y1n; y2n];
% 在此处计算Jacobi矩阵的值
ydn = double(subs(yd, [x1 x2], x_numeric));
end
%% 区间键入
function [area1, area2] = Area()
% 在此处定义根的搜索区间。若无法确定根的区间，请键入相同的任意两个数字。
area1 = [-10; -10];   % 变量的下限
area2 = [10; 10];     % 变量的上限
end
%% 牛顿法迭代流程
function [root, i, choose, wList] = Newton(max_step, max_f_step, max_x_step)
% 导入方程
[y1, y2, yd] = f();
% 计数
i = 1;
% 记录各步的root
wList = zeros(2, 1);
choose = 0;
% 取出区间端点
[area1, area2] = Area();
% w = [mean(area1); mean(area2)];  % 给定初始的root估计, 按列排布
w = [1.2; 0.8];
while i <= max_step
    wList(:, end+1) = w;
    % 计算函数值与导数值
    [yn, ydn] = f_calc(w', y1, y2, yd);
    % 函数/导数合理性检查
    if abs(max(yn)) == Inf || abs(max(max(ydn))) == Inf
        fprintf('警告! 迭代无法继续进行, 因为当前解下, 函数值或导数值趋于 +/- Inf. \n');
        disp('使用上一次(或初值)求出的根作为输出, 算法终止!');
        disp('请考虑使用非函数或导数的方法进行求解.');
        if i > 1
            root = wList(:, i-1);
        else
            root = w;
        end
        choose = 98;
        break
    end
    % 执行迭代
    wnew = w - ydn\yn;
    % 区间检查
    area1_check = wnew - area1;
    area2_check = area2 - wnew;
    if min(area1_check) < 0 || min(area2_check) < 0
        fprintf('警告！当前的根超出限定范围, 使用上一次(或初值)求出的根作为输出, 算法终止! \n');
        disp('请检查您规定的区间范围是否合理.');
        root = w;
        choose = 99;
        break
    end
    % 迭代停机检查
    if (norm(wnew-w) < max_x_step)
        root = wnew;
        choose = 1;
        break
    elseif (norm(f_calc(wnew', y1, y2, yd)) < max_f_step)
        root = wnew;
        choose = 2;
        break
    end
    w = wnew;
    i = i + 1;
end
end
%% 结果展示
function ShowResult(root, step, choose, wList)
disp('您提交的 使用牛顿法计算多元方程组的根 任务已经接近完成...');
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
    disp('迭代过程出现错误, 错误信息: 迭代过程中产生的函数值或导数值趋于无穷.');
end
fprintf('迭代求取方程的根为: \n');
disp(root);
[y1, y2, yd] = f();
fprintf('在此根下的函数值为: \n');
disp(f_calc(root', y1, y2, yd));
wList = wList(:, 2:end);
disp('打印 {步 - 误差} 与 {步 - log(误差)} 曲线...');
plot_x = 1:step;
plot_y = sqrt(sum((wList - root).^2));   % 相当于对各列取2-范数
plot_y2 = log(plot_y);
figure
subplot(1, 2, 1)
plot(plot_x, plot_y, '-*', 'LineWidth', 2, 'MarkerSize', 5);
set(gca,'FontSize', 15, 'Fontname', 'Times New Roman');
xlabel('Iteration Step', 'Fontname', 'Times New Roman','FontSize', 15);
ylabel('Error of Solution', 'Fontname', 'Times New Roman','FontSize', 15);
title('Iteration Result: Step - Error', 'Fontname', 'Times New Roman','FontSize', 18);
subplot(1, 2, 2)
plot(plot_x, plot_y2, '-*', 'Color', [1, 0.5, 0], 'LineWidth', 2, 'MarkerSize', 5);
set(gca,'FontSize', 15, 'Fontname', 'Times New Roman');
xlabel('Iteration Step', 'Fontname', 'Times New Roman','FontSize', 15);
ylabel('log(Error of Solution)', 'Fontname', 'Times New Roman','FontSize', 15);
title('Iteration Result: Step - log(Error)', 'Fontname', 'Times New Roman','FontSize', 18);
end