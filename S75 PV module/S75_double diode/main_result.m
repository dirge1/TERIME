clear;clc;close all
% 定义基础目录为当前目录
baseDir = pwd;

% 定义文件夹名称列表
folders = {'S_双二极管_200_25', 'S_双二极管_400_25', 'S_双二极管_600_25', 'S_双二极管_800_25', ...
    'S_双二极管_1000_25', 'S_双二极管_1000_20', 'S_双二极管_1000_40', 'S_双二极管_1000_60'};

% 遍历每一个文件夹
for w = 8
    folderPath = fullfile(baseDir, folders{w});
    
    % 获取文件夹中所有的 .mat 文件
    matFileNames{1} = 'result_RIME_improve_boundary2.mat';
    matFileNames{2} = 'result_RIME_differential_learning.mat';
    matFileNames{3} = 'result_DIWJAYA.mat';
    matFileNames{4} = 'result_Rao.mat';
    % 遍历每一个 .mat 文件
    for i = 1:length(matFileNames)
        fileName = fullfile(folderPath, matFileNames{i});
        
        % 从每个 .mat 文件中加载 result_RMSE
        data = load(fileName, 'result_RMSE','result_pa');
        result_RMSE=data.result_RMSE;
        result_pa=data.result_pa;
        [~, X] = min(result_RMSE); % Find the index of the minimum RMSE
        pab = result_pa{1}(4, X);  % Extract the parameters array
        % Extract each parameter and store in the matrix
        parameters{w}(i, 1) = pab{1}(1); % Iph
        parameters{w}(i, 2) = pab{1}(2)*1e6; % Io
        parameters{w}(i, 3) = pab{1}(3); % Rs
        parameters{w}(i, 4) = pab{1}(4); % Rsh
        parameters{w}(i, 5) = pab{1}(5); % n
        parameters{w}(i, 6) = pab{1}(6)*1e6; % Io2
        parameters{w}(i, 7) = pab{1}(7); % n2
        
        result(i,1)=min(result_RMSE)*100;
        result(i,2)=mean(result_RMSE)*100;
        result(i,3)=max(result_RMSE)*100;
        result(i,4)=std(result_RMSE);
    end
end
