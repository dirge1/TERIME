clear; clc; close all;

% List of result files
files = {
    'result_RIME_improve_boundary2.mat', 'result_RIME.mat', 'result_CDRIME.mat', ...
    'result_RIME_cross.mat', 'result_SRIME.mat', 'result_RIME_differential_learning.mat',...
    'result_DIWJAYA.mat', 'result_DO.mat', 'result_NGO.mat', 'result_RAO.mat'
};

% Initialize matrix to hold the parameters from each file
num_files = length(files);
parameters = zeros(num_files, 5); % 5 parameters per file

% Loop through each file
for i = 1:num_files
    % Load the file
    load(files{i});

    % Assuming 'result_pa' and 'result_RMSE' are loaded with the file
    [~, X] = min(result_RMSE); % Find the index of the minimum RMSE
    pab = result_pa{1}(4, X);  % Extract the parameters array

    % Extract each parameter and store in the matrix
    parameters(i, 1) = pab{1}(1); % Iph
    parameters(i, 2) = pab{1}(2)*1e6; % Io
    parameters(i, 3) = pab{1}(3); % Rs
    parameters(i, 4) = pab{1}(4); % Rsh
    parameters(i, 5) = pab{1}(5); % n
    parameters(i, 6) = pab{1}(6)*1e6; % Io2
    parameters(i, 7) = pab{1}(7); % n2
    
    result(i,1)=obj_min*1e3;
    result(i,2)=obj_mean*1e3;
    result(i,3)=obj_max*1e3;
    result(i,4)=std(result_RMSE);
    
end

