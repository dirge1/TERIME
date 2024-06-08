clear; clc; close all;
figure
% 数据
y = [1 1 1 1.375; 1 3 3.25 2.375; 1 3.375 3 2.25];
bar(y); % 增加条形边框宽度

% 设置图形尺寸和字体
fig_size = [10 5 20 15];
set(gcf, 'unit', 'centimeters', 'position', fig_size);
font_name = 'Times New Roman';
font_size = 20;

% 坐标轴设置
ax = gca;
set(ax, 'XTickLabel', {'Min', 'Mean', 'Max'}, 'fontsize', font_size, 'fontname', font_name);
xlabel('Parameter extraction on DDM for S75', 'fontsize', font_size);
ylabel('Average ranking', 'fontsize', font_size);

% 网格线和背景
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';
ax.GridColor = [0.8, 0.8, 0.8]; % 设定网格线颜色
ax.GridAlpha = 0.5; % 网格线透明度

% 图例
h = legend('TERIME', 'MRIME', 'DIWJAYA', 'CLRao-1', 'Location', 'northwest');
%% 
figure
% 数据
y=[1.37500000000000,1,2.62500000000000,1.62500000000000;1.25000000000000,1.75000000000000,4,2.25000000000000;1,2.37500000000000,3.75000000000000,2.12500000000000];
bar(y); % 增加条形边框宽度

% 设置图形尺寸和字体
fig_size = [10 5 20 15];
set(gcf, 'unit', 'centimeters', 'position', fig_size);
font_name = 'Times New Roman';
font_size = 20;

% 坐标轴设置
ax = gca;
set(ax, 'XTickLabel', {'Min', 'Mean', 'Max'}, 'fontsize', font_size, 'fontname', font_name);
xlabel('Parameter extraction on TDM for S75', 'fontsize', font_size);
ylabel('Average ranking', 'fontsize', font_size);

% 网格线和背景
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';
ax.GridColor = [0.8, 0.8, 0.8]; % 设定网格线颜色
ax.GridAlpha = 0.5; % 网格线透明度

% 图例
h = legend('TERIME', 'MRIME', 'DIWJAYA', 'CLRao-1', 'Location', 'northwest');
