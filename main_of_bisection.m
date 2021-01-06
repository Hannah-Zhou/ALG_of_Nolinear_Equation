%% ��׼�㷨1�����ַ�
% 
clear all;
clc;
%% ����ͣ��׼��
% ��������������
max_step = 200;
% ��������������޶�
max_f_step = 1e-16;
% ����������������޶�
max_x_step = 1e-14;
%% ������
[root, step, choose, wList] = Dichotomy(max_step, max_f_step, max_x_step);
%% ���չʾ
ShowResult(root, step, choose, wList);
%% ��Ԫ��������
function y = f(x)
% �ڴ˴����庯����Ҫ��Ϊ����Ԫ��
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
%% �������
function [a, b] = Area()
% �ڴ˴���������������䡣���޷�ȷ���������䣬�������ͬ�������������֡�
a = -0.4;
b = 0.71;
end
%% ���ַ���������
function [root, i, choose, wList] = Dichotomy(max_step, max_f_step, max_x_step)
% ����
i = 1;
% ��¼������root
wList = [];
choose = 0;
% ȡ������˵�
[a, b] = Area();
% �޸������������̶�
max_x_step = max_x_step + max([abs(a), abs(b), 1])*2^(-53);
% ��������
while i <= max_step
    w = (a + b)/2;
    wList(end+1) = w;
    % ����ͣ�����
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
%% ���չʾ
function ShowResult(root, step, choose, wList)
disp('���ύ�� ʹ�ö��ַ����㵥Ԫ���̵ĸ� �����Ѿ��ӽ����...');
disp('�㷨ִ�����, ��ӡ���...');
if choose == 0
    fprintf('������뺯��, �㷨�ڵ�%d����ɵ���, ʹ�õ�ͣ��׼��Ϊ�ﵽ������ָ����������. \n', step);
elseif choose == 1
    fprintf('������뺯��, �㷨�ڵ�%d����ɵ���, ʹ�õ�ͣ��׼��Ϊ�����С. \n', step);
else
    fprintf('������뺯��, �㷨�ڵ�%d����ɵ���, ʹ�õ�ͣ��׼��Ϊ����ֵ�������������޶�Ҫ��. \n', step);
end
fprintf('������ȡ���̵ĸ�Ϊ: %f. \n', root);
fprintf('�ڴ˸��µĺ���ֵΪ: %f. \n', f(root));
disp('��ӡ {�� - ���} �� {�� - log(���)} ����...');
plot_x = 1:step;
plot_y = abs(wList - root);
plot_x2 = log(plot_x);
plot_y2 = log(plot_y);
% �������
% �޳�-Inf���Inf��
CHfind_Index = find(plot_y2 ~= Inf & plot_y2 ~= -Inf);
linear_fit = polyfit(plot_x(CHfind_Index), plot_y2(CHfind_Index), 1);
plot_y2_fit = polyval(linear_fit, plot_x(CHfind_Index));
k = linear_fit(1,1);
disp_k = strcat('Linear k = ', num2str(k));
disp('���ܵ����������׿̻���');
fprintf('x - log(error)��������ֱ��б��k = %f. \n', k);
% �߽����������
plot_y3 = log(-plot_y2);
CHfind_Index_high = find(plot_y3 ~= Inf & plot_y3 ~= -Inf);
linear_fit_high = polyfit(plot_x(CHfind_Index_high), plot_y3(CHfind_Index_high), 1);
plot_y3_fit = polyval(linear_fit_high, plot_x(CHfind_Index_high));
k_high = linear_fit_high(1,1);
disp_k_high = strcat('Linear k = ', num2str(k_high));
disp('���ܵķ����������׿̻���');
fprintf('x - log(-log(error))��������ֱ��б��k = %f. \n', k_high);
% ͼ����ʾ
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





