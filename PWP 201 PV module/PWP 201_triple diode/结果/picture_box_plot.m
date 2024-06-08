clear;clc;close all
% Define the file names
files = {'result_RIME_improve_boundary2.mat', 'result_RIME_differential_learning.mat', 'result_DIWJAYA.mat', 'result_RAO.mat'};
labels = {'TERIME', 'MRIME', 'DIWJAYA', 'CLRao-1'};

% Initialize the data container
data = [];

% Loop through each file and load the result_RMSE
for i = 1:length(files)
    % Load the specific .mat file
    load(files{i});
    
    % Assuming 'result_RMSE' is the variable name holding the data
    if exist('result_RMSE', 'var')
        data = [data; result_RMSE];  % Append data
    else
        disp(['result_RMSE not found in ', files{i}]);
    end
end


% Create a box plot of the data
figure;
boxplot(data', 'Labels', labels);

title('Comparison of RMSE Across Algorithms');
ylabel('RMSE Values');

fosize=20;
set(gcf,'unit','centimeters','position',[10 5 25 15]);
xlabel('Iteration number','FontWeight','bold','fontsize',fosize,'fontname','Times New Roman');
ylabel('Mean RMSE','FontWeight','bold','fontsize',fosize,'fontname','Times New Roman');
set(gca, 'fontsize',fosize,'fontname','Times New Roman');