% To conduct graph theory analysis
% Xiao Chen
% 220411
%
% Redo the analysis of all Graph Theory results, especially betweenness
% Modified by Xiao Chen
% 241225
% 
% beta version: I did not reproduce the betweenness results. After redo the
% whole analysis pipeline, here is the beta version
% 250110
%
% Also do two sample t test
% 250507
%
% chenxiaophd@gmail.com

%% Initialization
clear; clc;
work_dir = '/mnt/Data3/RfMRILab/ChenX/rumination_network/analyses/redo_only_surface/graph_theory';
if ~exist(work_dir, 'dir'); mkdir(work_dir); end
data_dir = '/mnt/Data3/RfMRILab/ChenX/Suzhou_Rumination/Preprocessing';

GTA_metircs_rum_dir = ['/mnt/Data3/RfMRILab/ChenX/rumination_network/analyses', ...
                       '/redo_AN_FPN_DMN/GTA_metrics/rum'];
GTA_metrics_dis_dir = ['/mnt/Data3/RfMRILab/ChenX/rumination_network/analyses', ...
                       '/redo_AN_FPN_DMN/GTA_metrics/dis'];

% read in demographic data
demographic_info = readtable(['/mnt/Data3/RfMRILab/ChenX/rumination_network', ...
                              '/analyses/network_redo/sample_clean_imputed_v5.csv']);
sub_list = demographic_info.Serial_Number;
% group info, MDD or HC
Dx =  demographic_info.Dx;
Age = demographic_info.Age;
Sex = demographic_info.Sex;
Sex = Sex -1;
Edu = demographic_info.Edu;
HeadMotion = [demographic_info.head_motion_rum, demographic_info.head_motion_dis, ...
              demographic_info.head_motion_rest];

%% %%%%%%%%%% Analysis on rumination and distraction %%%%%%%%%%%%%%%%

%%%%%%%% Mixed effect Analysis %%%%%%%%
% Arrange data
% Rum
TempDir = [work_dir,'/Mixed_Effect/MDD_Rum'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = 1:length(find(Dx == 1))
    CurrentMatrixFile = [GTA_metircs_rum_dir, '/GTA_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end
TempDir = [work_dir,'/Mixed_Effect/HC_Rum'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = length(find(Dx == 1))+1:length(sub_list)
    CurrentMatrixFile = [GTA_metircs_rum_dir,'/GTA_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end

%Dis
TempDir = [work_dir,'/Mixed_Effect/MDD_Dis'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = 1:length(find(Dx == 1))
    CurrentMatrixFile = [GTA_metrics_dis_dir, '/GTA_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end
TempDir = [work_dir,'/Mixed_Effect/HC_Dis'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = length(find(Dx == 1))+1:length(sub_list)
    CurrentMatrixFile = [GTA_metrics_dis_dir, '/GTA_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end

%% do stats: mixed effect analysis

DependentDir{1,1} = [work_dir,'/Mixed_Effect/MDD_Rum'];
DependentDir{2,1} = [work_dir,'/Mixed_Effect/MDD_Dis'];
DependentDir{3,1} = [work_dir,'/Mixed_Effect/HC_Rum'];
DependentDir{4,1} = [work_dir,'/Mixed_Effect/HC_Dis'];

OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

% OtherCovariates = [];

mkdir([work_dir, '/Mixed_Effect/Stats']);
OutputName = [work_dir, '/Mixed_Effect/Stats/Mixed_Effect_onlycov_headmotion'];

y_MixedEffectsAnalysis_Image(DependentDir,OutputName,'',[],OtherCovariates, []);

%% paired t test
%MDD
DependentDir = [];OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/Mixed_Effect/MDD_Rum'];
DependentDir{2,1} = [work_dir,'/Mixed_Effect/MDD_Dis'];
OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
OutputDir = [work_dir, '/PairedT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir, '/PT_onlycov_headmotion_MDD'];

y_TTestPaired_Image(DependentDir,OutputName,'',[],OtherCovariates,[]);

% HC
DependentDir = [];OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/Mixed_Effect/HC_Rum'];
DependentDir{2,1} = [work_dir,'/Mixed_Effect/HC_Dis'];
OtherCovariates{1,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
OtherCovariates{2,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);
OutputDir = [work_dir, '/PairedT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir, '/PT_onlycov_headmotion_HC'];

y_TTestPaired_Image(DependentDir,OutputName,'',[],OtherCovariates,[]);

%% two sample t test

%rum
DependentDir = [];OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/Mixed_Effect/MDD_Rum'];
DependentDir{2,1} = [work_dir,'/Mixed_Effect/HC_Rum'];
OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
OtherCovariates{2,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
OutputDir = [work_dir, '/TwoT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir, '/TwoT_onlycov_headmotion_Rum'];

y_TTest2_Image(DependentDir,OutputName,'',[],OtherCovariates,[]);

%dis
DependentDir = [];OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/Mixed_Effect/MDD_Dis'];
DependentDir{2,1} = [work_dir,'/Mixed_Effect/HC_Dis'];
OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),2);
OtherCovariates{2,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);
OutputDir = [work_dir, '/TwoT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir, '/TwoT_onlycov_headmotion_Dis'];

y_TTest2_Image(DependentDir,OutputName,'',[],OtherCovariates,[]);


%% extract global metrics
MDD_Rum = []; MDD_Dis = [];
for j = 1:length(find(Dx == 1))
    load([DependentDir{1,1},'/GTA_sub-',sub_list{j},'.mat']);
    MDD_Rum(j,1) = Lambda_AUC(1); %change this
end
for j = 1:length(find(Dx == 1))
    load([DependentDir{2,1},'/GTA_sub-',sub_list{j},'.mat']);
    MDD_Dis(j,1) = Lambda_AUC(1); %change this
end

HC_Rum = []; HC_Dis = [];
count = 0;
for j = length(find(Dx == 1))+1:length(sub_list)
    count = count + 1;
    load([DependentDir{3,1},'/GTA_sub-',sub_list{j},'.mat']);
    HC_Rum(count,1) = Lambda_AUC(1); %change this
end
count = 0;
for j = length(find(Dx == 1))+1:length(sub_list)
    count = count + 1;
    load([DependentDir{4,1},'/GTA_sub-',sub_list{j},'.mat']);
    HC_Dis(count,1) = Lambda_AUC(1); %change this
end


% extract significant nodes' values
load([work_dir, '/Mixed_Effect/EigenvectorCentrality_NodeWeight.mat']);
index = find(NodeWeight);
temp = index(2);
MDD_Rum = []; MDD_Dis = [];
for j = 1:length(find(Dx == 1))
    load([DependentDir{1,1},'/GTA_sub-',sub_list{j},'.mat']);
    MDD_Rum(j,1) = EigenvectorCentrality_AUC(temp); 
end
for j = 1:length(find(Dx == 1))
    load([DependentDir{2,1},'/GTA_sub-',sub_list{j},'.mat']);
    MDD_Dis(j,1) = EigenvectorCentrality_AUC(temp); 
end

HC_Rum = []; HC_Dis = [];
count = 0;
for j = length(find(Dx == 1))+1:length(sub_list)
    count = count + 1;
    load([DependentDir{3,1},'/GTA_sub-',sub_list{j},'.mat']);
    HC_Rum(count,1) = EigenvectorCentrality_AUC(temp); 
end
count = 0;
for j = length(find(Dx == 1))+1:length(sub_list)
    count = count + 1;
    load([DependentDir{4,1},'/GTA_sub-',sub_list{j},'.mat']);
    HC_Dis(count,1) = EigenvectorCentrality_AUC(temp); 
end