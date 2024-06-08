clear; clc; close all;
figure
% 数据
y = [1 1 1 1; 1 3 2 4; 1 4 2 3];
bar(y); % 增加条形边框宽度

% 设置图形尺寸和字体
fig_size = [10 5 20 15];
set(gcf, 'unit', 'centimeters', 'position', fig_size);
font_name = 'Times New Roman';
font_size = 20;

% 坐标轴设置
ax = gca;
set(ax, 'XTickLabel', {'Min', 'Mean', 'Max'}, 'fontsize', font_size, 'fontname', font_name);
xlabel('Parameter extraction on DDM for RTC France', 'fontsize', font_size);
ylabel('Ranking', 'fontsize', font_size);

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
y=[2 1 4 3; 1 2 3 4; 1 5 3 6];
bar(y); % 增加条形边框宽度

% 设置图形尺寸和字体
fig_size = [10 5 20 15];
set(gcf, 'unit', 'centimeters', 'position', fig_size);
font_name = 'Times New Roman';
font_size = 20;

% 坐标轴设置
ax = gca;
set(ax, 'XTickLabel', {'Min', 'Mean', 'Max'}, 'fontsize', font_size, 'fontname', font_name);
xlabel('Parameter extraction on TDM for RTC France', 'fontsize', font_size);
ylabel('Ranking', 'fontsize', font_size);

% 网格线和背景
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';
ax.GridColor = [0.8, 0.8, 0.8]; % 设定网格线颜色
ax.GridAlpha = 0.5; % 网格线透明度

% 图例
h = legend('TERIME', 'MRIME', 'DIWJAYA', 'CLRao-1', 'Location', 'northwest');
