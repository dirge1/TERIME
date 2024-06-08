clear;clc;close all
% �������Ŀ¼Ϊ��ǰĿ¼
baseDir = pwd;

% �����ļ��������б�
folders = {'S_˫������_200_25', 'S_˫������_400_25', 'S_˫������_600_25', 'S_˫������_800_25', ...
    'S_˫������_1000_25', 'S_˫������_1000_20', 'S_˫������_1000_40', 'S_˫������_1000_60'};
labels = {'TERIME', 'MRIME', 'DIWJAYA', 'CLRao-1'};
% ��ʼ��һ���ṹ�����洢 RMSE ����
results = struct();
data = [];
% ����ÿһ���ļ���
for w = 8
    folderPath = fullfile(baseDir, folders{w});
    
    % ��ȡ�ļ��������е� .mat �ļ�
    matFileNames{1} = 'result_RIME_improve_boundary2.mat';
    matFileNames{2} = 'result_RIME_differential_learning.mat';
    matFileNames{3} = 'result_DIWJAYA.mat';
    matFileNames{4} = 'result_Rao.mat';
    % ����ÿһ�� .mat �ļ�
    for i = 1:length(matFileNames)
        fileName = fullfile(folderPath, matFileNames{i});
        
        % ��ÿ�� .mat �ļ��м��� result_RMSE
        data2 = load(fileName, 'result_RMSE','result_pa');
        result_RMSE=data2.result_RMSE;
        data = [data; result_RMSE];  % Append data
    end
end


% Create a box plot of the data
figure;
boxplot(data', 'Labels', labels);

title('Comparison of RMSE under 1000W/m^2 60��');
ylabel('RMSE Values');

fosize=20;
set(gcf,'unit','centimeters','position',[10 5 25 15]);
xlabel('Iteration number','FontWeight','bold','fontsize',fosize,'fontname','Times New Roman');
ylabel('Mean RMSE','FontWeight','bold','fontsize',fosize,'fontname','Times New Roman');
set(gca, 'fontsize',fosize,'fontname','Times New Roman');