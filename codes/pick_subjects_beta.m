% To pick the subjects to do the beta version of analyses
%
% Written by Xiao Chen
% 241230
%
% Revised: merging behavior ratings in the final demograohic info table
% Xiao Chen
% 250124
%
% chenxiaophd@gmail.com

%% initialization
clear; clc;

work_dir = '/mnt/Data3/RfMRILab/ChenX/rumination_network/analyses/network_redo';
if ~exist(work_dir, 'dir'); mkdir(work_dir); end
data_dir = '/mnt/Data3/RfMRILab/ChenX/Suzhou_Rumination/Preprocessing';

%% read in the demographic data and do cleaning
demographic_info = readtable([work_dir, '/demographic_info_v10.xlsx'], 'Sheet', 1);
% select some columns
selected_var = {'Serial_Number', 'Dx', 'Sex', 'Age', 'Edu', 'Handedness_label', 'Drug_naive', ...
                'HAMD_Sum', 'HAMA_Sum', 'Rumination', 'Reflection', 'Brooding', 'BDI'};
selected_data = demographic_info(:, selected_var);
% set all 999s as NaN
for i = 1:width(selected_data)
    if isnumeric(selected_data{:, i})
        selected_data{selected_data{:, i} == 999, i} = NaN;
        selected_data{selected_data{:, i} == 888, i} = NaN; %some variables NA to HCs
    end
end
% round all imputed data to the nearest integers
selected_var = {'HAMD_Sum', 'HAMA_Sum', 'Rumination', 'Reflection', 'Brooding', 'BDI'};
for i = 1:length(selected_var)
    columnName = selected_var{i};
    selected_data.(columnName) = round(selected_data.(columnName));
end

% read in the behavior data
behavior_data = readtable([work_dir, '/behaviral_results.csv']);
% set all 999s as NaN
for i = 1:width(behavior_data)
    if isnumeric(behavior_data{:, i})
        behavior_data{behavior_data{:, i} == 999, i} = NaN;
        behavior_data{behavior_data{:, i} == 888, i} = NaN; %some variables NA to HCs
    end
end

% do exclusion
% exclude some un-qualified subjects
% #1: if age, sex, and Edu has any NaNs
rowsToExclude = isnan(selected_data.Age) | isnan(selected_data.Sex) | isnan(selected_data.Edu);
selected_data = selected_data(~rowsToExclude, :);

% #2: exclude any left-handed
selected_data = selected_data(selected_data.Handedness_label ~= 2, :);

% #3: exclude two subjects with wrong phase coding direction
selected_data = selected_data(~ismember(selected_data.Serial_Number, {'Sub026', 'Sub051'}), :);

% #4: have all the resting state, rumination/distraction states data
sub_list_task = importdata([work_dir, '/sub_list_task_merged20220118.txt']);
sub_list_rest = importdata([work_dir, '/sub_list_rest_merged20220105.txt']);
commonSubjects = intersect(sub_list_task, sub_list_rest);
selected_data = selected_data(ismember(selected_data.Serial_Number, commonSubjects), :);

% #5: rule out excessive head motion, FD_Jenkinson > 0.2
% 1: Rum, 2: Dis, 3: Rest
HeadMotion = zeros(height(selected_data),3);
for i=1:height(selected_data)
    Temp=load([data_dir,'/Task/RealignParameter/sub-',selected_data.Serial_Number{i}, ...
                                '/FD_Jenkinson_sub-',selected_data.Serial_Number{i},'.txt']);
    HeadMotion(i,1)=mean(Temp);
     Temp=load([data_dir,'/Task/RealignParameter/sub-',selected_data.Serial_Number{i}, ...
                             '/S2_FD_Jenkinson_sub-',selected_data.Serial_Number{i},'.txt']);
    HeadMotion(i,2)=mean(Temp);
    Temp=load([data_dir,'/Rest/RealignParameter/sub-',selected_data.Serial_Number{i}, ...
                             '/FD_Jenkinson_sub-',selected_data.Serial_Number{i},'.txt']);
    HeadMotion(i,3)=mean(Temp);
end
matrixTable = array2table(HeadMotion, ...
            'VariableNames', {'head_motion_rum', 'head_motion_dis', 'head_motion_rest'});
selected_data = [selected_data, matrixTable];
rowsToExclude = selected_data.head_motion_rum > 0.2 | ...
                    selected_data.head_motion_dis > 0.2 | ...
                    selected_data.head_motion_rest > 0.2;
selected_data = selected_data(~rowsToExclude, :);

%% do exclusion to balance
% To balance the patient and HC samples
selected_data = selected_data(~ismember(selected_data.Serial_Number, ...
                                {'Sub058','Sub031', 'Sub139','Sub123','Sub116'}), :);
                            
% selected_data = selected_data(~ismember(selected_data.Serial_Number, ...
%                                 {'Sub116','Sub123', 'Sub139', 'Sub151','Sub139'}), :);

%% merging the behavior data with the demographic data
behavior_data.Properties.VariableNames{1} = 'Serial_Number';
selected_merged_data = join(selected_data, behavior_data, 'Keys', 'Serial_Number');

%% test the balance of Age, Sex ratio and edu
% Subset the data into two groups based on Dx
group1 = selected_data(selected_data.Dx == 1, :);  % Dx == 1
group2 = selected_data(selected_data.Dx == 2, :);  % Dx == 2

n_group1 = height(group1);  % Sample size of group 1
n_group2 = height(group2);  % Sample size of group 2
fprintf('Sample size for Dx == 1: %d\n', n_group1);
fprintf('Sample size for Dx == 2: %d\n\n', n_group2);

% 1. Independent t-test for Age
[h_age, p_age, ci_age, stats_age] = ttest2(group1.Age, group2.Age);
fprintf('Independent t-test for Age:\n');
fprintf('t-statistic = %.4f\n', stats_age.tstat);
fprintf('p-value = %.4f\n', p_age);
fprintf('Confidence interval = [%.4f, %.4f]\n\n', ci_age(1), ci_age(2));

% 2. Independent t-test for Edu
[h_edu, p_edu, ci_edu, stats_edu] = ttest2(group1.Edu, group2.Edu);
fprintf('Independent t-test for Edu:\n');
fprintf('t-statistic = %.4f\n', stats_edu.tstat);
fprintf('p-value = %.4f\n', p_edu);
fprintf('Confidence interval = [%.4f, %.4f]\n\n', ci_edu(1), ci_edu(2));

% 3. Chi-square test for Sex (contingency table)
% Create a contingency table for Sex by Dx
[~,chi2,p,~] = crosstab(selected_data.Dx, selected_data.Sex);
fprintf('Chi-square test for Sex distribution:\n');
fprintf('Chi-square statistic = %.4f\n', chi2);
fprintf('p-value = %.4f\n', p);


%% test whether the final sample is feasible
% See is there any remaining NaNs
nanColumns = any(ismissing(selected_data), 1);
columnsWithNaNs = selected_data.Properties.VariableNames(nanColumns);
if isempty(columnsWithNaNs)
    disp('No columns contain NaN values.');
else
    disp('Columns containing NaN values:');
    disp(columnsWithNaNs);
end

% to get the potential to-be-excluded participants
group2.benchmark = 0*group2.Age + 1*group2.Edu;
[~, idx] = max(group2.benchmark);  
subject_interested = group2(idx, :);

% to get the potential to-be-excluded participants
group1.benchmark = 0.5*group1.Age + 0.5*group1.Edu;
[~, idx] = min(group1.benchmark);  
subject_interested = group1(idx, :);

%% write out the final analyzed sample
% v2: also included drug_naive
% v3: exclude Sub107 to balance the age
% v4: change the excluded subject list to balance better
% v5: merging the behavior data to the demographic data
writetable(selected_merged_data, [work_dir, '/sample_cleaned_merged_v5.csv']);