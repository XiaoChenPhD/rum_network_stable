% explore the correlation between the graph theory metrics and the scales,
% as well as the behavioral ratings
%
% Xiao Chen
% 250303
% chenxiaophd@gmail.com

%% initialization

clear;clc;

work_dir = 'xxxxxxxx';
stats_dir = 'xxxxxxxx';

demographic_info = readtable('xxxxxxxxx/sample_clean_imputed_v5.csv');
load([stats_dir, '/mixed_effect_interaction_stats_corrected.mat']);
Dx = demographic_info.Dx;
sub_list = demographic_info.Serial_Number;
                     
%% get the interested GTA metrics and scales

scale_set = {'HAMD_Sum'};
scale_matrix = demographic_info{:, scale_set};
subj_idx = find(demographic_info.Dx == 1); %only explore correlation in the MDD patients


% visually check the significant results
sig_results = stats_output.EigenvectorCentrality_AUC;
idx_current = 295;%use the atlas index instead
label_current = sig_results{2,5};

% MDD
sub_list_MDD = sub_list(Dx == 1);
rum_metric_MDD = zeros(length(sub_list_MDD),1);
dis_metric_MDD = zeros(length(sub_list_MDD),1);
for i = 1:length(sub_list_MDD)
    % rum
    load('xxxxx/GTA_',sub_list_MDD{i}, '.mat');
    rum_metric_MDD(i) = EigenvectorCentrality_AUC(idx_current);
    
    % dis
    load('xxxxx/GTA_',sub_list_MDD{i}, '.mat');
    dis_metric_MDD(i) = EigenvectorCentrality_AUC(idx_current);
end


% HC
sub_list_HC = sub_list(Dx == 2);
rum_metric_HC = zeros(length(sub_list_HC),1);
dis_metric_HC = zeros(length(sub_list_HC),1);
for i = 1:length(sub_list_HC)
    % rum
    load(['xxxxx/GTA_',sub_list_HC{i}, '.mat']);
    rum_metric_HC(i) = EigenvectorCentrality_AUC(idx_current);
    
    % dis
    load(['xxxxxx/GTA_',sub_list_HC{i}, '.mat']);
    dis_metric_HC(i) = EigenvectorCentrality_AUC(idx_current);
end

%% do corr analysis

delta_MDD = rum_metric_MDD - dis_metric_MDD;
delta_relative_MDD = (rum_metric_MDD - dis_metric_MDD)./dis_metric_MDD;

for j = 1:length(scale_set)
    % each field consists of a n x 4 matrix
    % sig index, sig label, r values, p values
    r_rum.(scale_set{j}) = cell(length(idx_current),4);
    r_delta.(scale_set{j}) = cell(length(idx_current),4);
    r_delta_relative.(scale_set{j}) = cell(length(idx_current),4);
end


for j = 1:size(scale_matrix,2)
    y = scale_matrix(subj_idx, j);
    for i = 1:length(idx_current)
        r_rum.(scale_set{j}){i,1} = idx_current;
        r_rum.(scale_set{j}){i,2} = label_current;
        x_rum = rum_metric_MDD(:,i);
        r_delta.(scale_set{j}){i,1} = idx_current;
        r_delta.(scale_set{j}){i,2} = label_current;
        x_delta = delta_MDD(:,i);
        r_delta_relative.(scale_set{j}){i,1} = idx_current;
        r_delta_relative.(scale_set{j}){i,2} = label_current;
        x_delta_relative = delta_relative_MDD(:,i);
        [r_rum.(scale_set{j}){i,3}, r_rum.(scale_set{j}){i,4}] = corr(x_rum, y);
        [r_delta.(scale_set{j}){i,3}, r_delta.(scale_set{j}){i,4}] = corr(x_delta, y);
        [r_delta_relative.(scale_set{j}){i,3}, r_delta_relative.(scale_set{j}){i,4}] = corr(x_delta_relative, y);
    end
end

output_name = 'GTA_corr';
save([work_dir, '/', output_name], ...
    'r_rum', 'r_delta', 'r_delta_relative');

%% extract numbers for plot

findings_scale = {'HAMD_Sum'};
scale_vec = demographic_info{:, findings_scale{1}}(Dx == 1);
matrix_output = [delta_MDD, scale_vec];

output_name = [label_current, '_', findings_scale{1}];
writematrix(matrix_output, [work_dir, '/', output_name, '.csv']);