clear;clc;close all
% 定义基础目录为当前目录
baseDir = pwd;

% 定义文件夹名称列表
folders = {'S_三二极管_200_25', 'S_三二极管_400_25', 'S_三二极管_600_25', 'S_三二极管_800_25', ...
    'S_三二极管_1000_25', 'S_三二极管_1000_20', 'S_三二极管_1000_40', 'S_三二极管_1000_60'};
labels = {'TERIME', 'MRIME', 'DIWJAYA', 'CLRao-1'};
% 初始化一个结构体来存储 RMSE 数据
results = struct();
data = [];
% 遍历每一个文件夹
for w = 5
    folderPath = fullfile(baseDir, folders{w});
    
    % 获取文件夹中所有的 .mat 文件
    matFiles = dir(fullfile(folderPath, '*.mat'));
    matFileNames2 = {matFiles.name};
    matFileNames{1} = matFileNames2{4};
    matFileNames{2} = matFileNames2{3};
    matFileNames{3} = matFileNames2{1};
    matFileNames{4} = matFileNames2{2};
    % 遍历每一个 .mat 文件
    for i = 1:length(matFileNames)
        fileName = fullfile(folderPath, matFileNames{i});
        
        % 从每个 .mat 文件中加载 result_RMSE
        data2 = load(fileName, 'result_RMSE','result_pa');
        result_RMSE=data2.result_RMSE;
        data = [data; result_RMSE];  % Append data
    end
end


% Create a box plot of the data
figure;
boxplot(data', 'Labels', labels);

title('Comparison of RMSE under 600W/m^2 25℃');
ylabel('RMSE Values');

fosize=20;
set(gcf,'unit','centimeters','position',[10 5 25 15]);
xlabel('Iteration number','FontWeight','bold','fontsize',fosize,'fontname','Times New Roman');
ylabel('Mean RMSE','FontWeight','bold','fontsize',fosize,'fontname','Times New Roman');
set(gca, 'fontsize',fosize,'fontname','Times New Roman');