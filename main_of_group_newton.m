%% ����ʵ�ַ�������ʽ��ţ�ٵ�������
% ������ʵ�ֻ���ţ�ٵ�����.
% Write by Zhou Ruihan
clear all;
clc;
%% ����ͣ��׼��
% ��������������
max_step = 1000;
% ��������������޶�
max_f_step = 1e-16;
% ����������������޶�
max_x_step = 1e-16;
%% ������
[root, step, choose, wList] = Newton(max_step, max_f_step, max_x_step);
%% ���չʾ
ShowResult(root, step, choose, wList);

%% ��Ԫ��������
function [y1, y2, yd] = f()
% �ڴ˴����庯����Ҫ��Ϊ���Ԫ��
syms x1 x2 x3
y1 = exp(x1^2 + x2^2) - exp(2);
y2 = exp(x1^2 - x2^2) - 1;
% ����Jacobi����
yd = jacobian([y1; y2], [x1 x2]);
end
%% ���㺯��ֵ��Jacobi�����ֵ
function [yn, ydn] = f_calc(x_numeric, y1, y2, yd)
% �ڴ˴����㺯��ֵ
syms x1 x2 x3
y1n = double(subs(y1, [x1 x2], x_numeric));
y2n = double(subs(y2, [x1 x2], x_numeric));
yn = [y1n; y2n];
% �ڴ˴�����Jacobi�����ֵ
ydn = double(subs(yd, [x1 x2], x_numeric));
end
%% �������
function [area1, area2] = Area()
% �ڴ˴���������������䡣���޷�ȷ���������䣬�������ͬ�������������֡�
area1 = [-10; -10];   % ����������
area2 = [10; 10];     % ����������
end
%% ţ�ٷ���������
function [root, i, choose, wList] = Newton(max_step, max_f_step, max_x_step)
% ���뷽��
[y1, y2, yd] = f();
% ����
i = 1;
% ��¼������root
wList = zeros(2, 1);
choose = 0;
% ȡ������˵�
[area1, area2] = Area();
% w = [mean(area1); mean(area2)];  % ������ʼ��root����, �����Ų�
w = [1.2; 0.8];
while i <= max_step
    wList(:, end+1) = w;
    % ���㺯��ֵ�뵼��ֵ
    [yn, ydn] = f_calc(w', y1, y2, yd);
    % ����/���������Լ��
    if abs(max(yn)) == Inf || abs(max(max(ydn))) == Inf
        fprintf('����! �����޷���������, ��Ϊ��ǰ����, ����ֵ����ֵ���� +/- Inf. \n');
        disp('ʹ����һ��(���ֵ)����ĸ���Ϊ���, �㷨��ֹ!');
        disp('�뿼��ʹ�÷Ǻ��������ķ����������.');
        if i > 1
            root = wList(:, i-1);
        else
            root = w;
        end
        choose = 98;
        break
    end
    % ִ�е���
    wnew = w - ydn\yn;
    % ������
    area1_check = wnew - area1;
    area2_check = area2 - wnew;
    if min(area1_check) < 0 || min(area2_check) < 0
        fprintf('���棡��ǰ�ĸ������޶���Χ, ʹ����һ��(���ֵ)����ĸ���Ϊ���, �㷨��ֹ! \n');
        disp('�������涨�����䷶Χ�Ƿ����.');
        root = w;
        choose = 99;
        break
    end
    % ����ͣ�����
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
%% ���չʾ
function ShowResult(root, step, choose, wList)
disp('���ύ�� ʹ��ţ�ٷ������Ԫ������ĸ� �����Ѿ��ӽ����...');
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
    disp('�������̳��ִ���, ������Ϣ: ���������в����ĺ���ֵ����ֵ��������.');
end
fprintf('������ȡ���̵ĸ�Ϊ: \n');
disp(root);
[y1, y2, yd] = f();
fprintf('�ڴ˸��µĺ���ֵΪ: \n');
disp(f_calc(root', y1, y2, yd));
wList = wList(:, 2:end);
disp('��ӡ {�� - ���} �� {�� - log(���)} ����...');
plot_x = 1:step;
plot_y = sqrt(sum((wList - root).^2));   % �൱�ڶԸ���ȡ2-����
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