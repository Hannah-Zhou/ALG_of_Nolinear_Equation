/%% ��׼�㷨4����ʹ�õ�����Ϣ�����ʹ��һ��C0��Ϣ�Ļ��� - setff�����㷨
% Write by Zhou Ruihan
% Date: 2020 / 10 / 23
% Version: 0.1 Trial V.
% BASIC VERSION: (1)The int - steff algorithm without derivative information is proposed.
clear all;
clc;
%% ����ͣ��׼��
% ��������������
max_step = 200;
% ��������������޶�
max_f_step = 1e-16;
% ����������������޶�
max_x_step = 1e-14;
% ����������
start_point = 0.041;
%% ������
[root, step, choose, wList] = IntSteff(max_step, max_f_step, max_x_step, start_point);
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
% �ڴ˴���������������䣬���˱�������෴��
a = -0.4;
b = 0.71;
end
%% Steffensen�㷨��������
function [root, i, choose, wList] = IntSteff(max_step, max_f_step, max_x_step, start_point)
% ����
i = 1;
% ȡ������˵�
[a, b] = Area();
% w = start_point;
w = (a + b) / 2;
% ��¼������root
wList = [w];
choose = 0;
% �޸������������̶�
max_x_step = max_x_step + max([abs(a), abs(b), 1])*2^(-53);
% ��ʼ����
while i <= max_step
    % ����/���������Լ�� 
    if abs(f(w)) == Inf
        fprintf('����! �����޷���������, ��Ϊ��ǰ��: root = %f��, ����ֵ�������ֵ�������� +/- Inf. \n', w);
        disp('ʹ����һ��(���ֵ)����ĸ���Ϊ���, �㷨��ֹ!');
        disp('�뿼��ʹ�÷Ǻ����ķ����������.');
        root = w;
        choose = 98;
        break
    end
    % ִ�е���
    alpha = 0.5;   % ����alpha��ֵ
    z = f(w);
    d = w - alpha * z^2 / (f(w + z) - z);
    wnew = w + 3*(d - w)*z / (2*z - f(d));
    if abs(wnew) == Inf
        fprintf('���棡��ǰ�ĵ�������ɢ: %f, ʹ����һ��(���ֵ)����ĸ���Ϊ���, �㷨��ֹ! \n', w);
        disp('�������涨�����䷶Χ�ͳ�ʼ�������Ƿ����.');
        root = w;
        choose = 97;
        break
    end
    % ������
    if wnew > b || wnew < a
        fprintf('���棡��ǰ�ĸ������޶���Χ: %f, ʹ����һ��(���ֵ)����ĸ���Ϊ���, �㷨��ֹ! \n', w);
        disp('�������涨�����䷶Χ�Ƿ����.');
        root = w;
        choose = 99;
        break
    end
    % ����ͣ�����
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
%% ���չʾ
function ShowResult(root, step, choose, wList)
disp('���ύ�� ʹ�û���-Steff�����㷨���㵥Ԫ���̵ĸ� �����Ѿ��ӽ����...');
disp('�㷨ִ�����, ��ӡ���...');
if choose == 0
    fprintf('������뺯��, �㷨�ڵ�%d����ɵ���, ʹ�õ�ͣ��׼��Ϊ�ﵽ������ָ����������. \n', step);
elseif choose == 1
    fprintf('������뺯��, �㷨�ڵ�%d����ɵ���, ʹ�õ�ͣ��׼��Ϊ���ڵĵ���������С, ��������. \n', step);
elseif choose == 2
    fprintf('������뺯��, �㷨�ڵ�%d����ɵ���, ʹ�õ�ͣ��׼��Ϊ����ֵ�������������޶�Ҫ��. \n', step);
elseif choose == 99
    disp('�������̳��ִ���, ������Ϣ: ���������Ľⳬ�����޶������䷶Χ.');
elseif choose == 98
    disp('�������̳��ִ���, ������Ϣ: ���������в����ĺ���ֵ�����ֵ������������.');
elseif choose == 97
    disp('�������̳��ִ���, ������Ϣ: ���������г����˵�������ɢ���������ʹ�õ��������޷�����.')
end
fprintf('������ȡ���̵ĸ�Ϊ: %f. \n', root);
fprintf('�ڴ˸��µĺ���ֵΪ: %f. \n', f(root));
disp('��ӡ {�� - ���} �� {�� - log(���)} ����...');
plot_x = 1:step;
plot_y = abs(wList(1:end) - root);
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