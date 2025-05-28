% To conduct Network statistic analysis for Suzhou Rumination Project
% Written by ChenXiao
% 20220302
% 
% Redo the whole brain (Schaefer 400 17 networks + subcortical 54) analyses
% Using NBS as multiple comparsion correction strategy
% Modified by Xiao Chen
% 20241213
%
% beta version
% 250102
%
% add a one sample t 
% 250318
%
% chenxiaophd@gmail.com

%% initialization
clear; clc;
work_dir = '/mnt/Data3/RfMRILab/ChenX/rumination_network/analyses/redo_only_surface/do_comparison';
if ~exist(work_dir, 'dir'); mkdir(work_dir); end
network_dir = '/mnt/Data3/RfMRILab/ChenX/rumination_network/analyses/redo_only_surface/network_matrix';
data_dir = '/mnt/Data3/RfMRILab/ChenX/Suzhou_Rumination/Preprocessing';

%% initilize the location of network matrixes
network_rum_dir = [network_dir, '/network_level/rum/'];
nodal_rum_dir = [network_dir, '/nodal_level/rum/'];
network_dis_dir = [network_dir, '/network_level/dis/'];
nodal_dis_dir = [network_dir, '/nodal_level/dis/'];
network_rest_dir = [network_dir, '/network_level/rest/'];
nodal_rest_dir = [network_dir, '/nodal_level/rest/'];
% network_sad_dir = [network_dir, '/Network_Level/NetworkMatrix_Sad/'];
% nodal_sad_dir = [network_dir, '/Nodal_Level/NetworkMatrix_Sad/'];

%% read in demographic data
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

%%  subgroup comparisons
% % include only drugnaive subjects
% WantedSubMatrix = ones(length(sub_list),1);
% Drug = demographic_info.Drug_naive;
% WantedSubMatrix(find(Drug == 0)) = 0;
% WantedSubMatrix(find(Drug == 999)) = 0;
% WantedSubIndex = find(WantedSubMatrix);
% sub_list = sub_list(WantedSubIndex);
% Dx = Dx(WantedSubIndex);
% Age = Age(WantedSubIndex);
% Sex = Sex(WantedSubIndex);
% Edu = Edu(WantedSubIndex);

% % include only recurrent subjects
% WantedSubMatrix = ones(length(sub_list),1);
% Drug = demographic_info.Drug_naive;
% WantedSubMatrix(find(Drug == 1)) = 0;
% WantedSubMatrix(find(Drug == 999)) = 0;
% WantedSubIndex = find(WantedSubMatrix);
% sub_list = sub_list(WantedSubIndex);
% Dx = Dx(WantedSubIndex);
% Age = Age(WantedSubIndex);
% Sex = Sex(WantedSubIndex);
% Edu = Edu(WantedSubIndex);

%% create a log file

currentDate = datestr(now, 'yyyy-mm-dd');
log_file = sprintf([work_dir, '/log_%s.txt'], currentDate);
logFileID = fopen(log_file, 'a');

%% Arrange data
% Network level
% Rum
TempDir = [work_dir,'/NetworkLevel/Mixed_Effect/MDD_Rum'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = 1:length(find(Dx == 1))
    CurrentMatrixFile = [network_rum_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end

TempDir = [work_dir,'/NetworkLevel/Mixed_Effect/HC_Rum'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = length(find(Dx == 1))+1:length(sub_list)
    CurrentMatrixFile = [network_rum_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end

%Dis
TempDir = [work_dir,'/NetworkLevel/Mixed_Effect/MDD_Dis'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = 1:length(find(Dx == 1))
    CurrentMatrixFile = [network_dis_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end

TempDir = [work_dir,'/NetworkLevel/Mixed_Effect/HC_Dis'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = length(find(Dx == 1))+1:length(sub_list)
    CurrentMatrixFile = [network_dis_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end

%rest
TempDir = [work_dir,'/NetworkLevel/Mixed_Effect/MDD_Rest'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = 1:length(find(Dx == 1))
    CurrentMatrixFile = [network_rest_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end

TempDir = [work_dir,'/NetworkLevel/Mixed_Effect/HC_Rest'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = length(find(Dx == 1))+1:length(sub_list)
    CurrentMatrixFile = [network_rest_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end

% Arrange data
% Rum
TempDir = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Rum'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = 1:length(find(Dx == 1))
    CurrentMatrixFile = [nodal_rum_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end
TempDir = [work_dir,'/NodalLevel/Mixed_Effect/HC_Rum'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = length(find(Dx == 1))+1:length(sub_list)
    CurrentMatrixFile = [nodal_rum_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end

%Dis
TempDir = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Dis'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = 1:length(find(Dx == 1))
    CurrentMatrixFile = [nodal_dis_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end
TempDir = [work_dir,'/NodalLevel/Mixed_Effect/HC_Dis'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = length(find(Dx == 1))+1:length(sub_list)
    CurrentMatrixFile = [nodal_dis_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end

%Rest
TempDir = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Rest'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = 1:length(find(Dx == 1))
    CurrentMatrixFile = [nodal_rest_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end
TempDir = [work_dir,'/NodalLevel/Mixed_Effect/HC_Rest'];
if ~exist(TempDir, 'dir')
    mkdir(TempDir); 
end
for i = length(find(Dx == 1))+1:length(sub_list)
    CurrentMatrixFile = [nodal_rest_dir,'NetworkMatrix_',sub_list{i},'.mat'];
    copyfile(CurrentMatrixFile,TempDir);
end

%% do stats
% network level
DependentDir{1,1} = [work_dir,'/NetworkLevel/Mixed_Effect/MDD_Rum'];
DependentDir{2,1} = [work_dir,'/NetworkLevel/Mixed_Effect/MDD_Dis'];
DependentDir{3,1} = [work_dir,'/NetworkLevel/Mixed_Effect/HC_Rum'];
DependentDir{4,1} = [work_dir,'/NetworkLevel/Mixed_Effect/HC_Dis'];

% also include age, sex, edu as covariates
% OtherCovariates{1,1} = [Age(1:length(find(Dx == 1)),1), Sex(1:length(find(Dx == 1)),1), Edu(1:length(find(Dx == 1)),1), HeadMotion(1:length(find(Dx == 1)),1)];
% OtherCovariates{2,1} = [Age(1:length(find(Dx == 1)),1), Sex(1:length(find(Dx == 1)),1), Edu(1:length(find(Dx == 1)),1), HeadMotion(1:length(find(Dx == 1)),2)];
% OtherCovariates{3,1} = [Age(length(find(Dx == 1))+1:length(sub_list),1), Sex(length(find(Dx == 1))+1:length(sub_list),1), Edu(length(find(Dx == 1))+1:length(sub_list),1), HeadMotion(length(find(Dx == 1))+1:length(sub_list),1)];
% OtherCovariates{4,1} = [Age(length(find(Dx == 1))+1:length(sub_list),1), Sex(length(find(Dx == 1))+1:length(sub_list),1), Edu(length(find(Dx == 1))+1:length(sub_list),1), HeadMotion(length(find(Dx == 1))+1:length(sub_list),2)];

% only include head motion as covariate
OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

% no covariates
% OtherCovariates = [];

% % PALM settings
% PALMSettings.nPerm = 5000;
% PALMSettings.ClusterInference=0;
% PALMSettings.ClusterFormingThreshold=2.3;
% PALMSettings.TFCE=0;
% PALMSettings.FDR=0;
% PALMSettings.TwoTailed=1;
% PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
% PALMSettings.SavePermutations = 1;

mkdir([work_dir, '/NetworkLevel/Mixed_Effect/Stats']);
OutputName = [work_dir, '/NetworkLevel/Mixed_Effect/Stats/Mixed_Effect_FDR'];


try
    y_MixedEffectsAnalysis_Image(DependentDir,OutputName,'',[],OtherCovariates, []);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end

%%
%%%%%%%%% Nodal Level %%%%%%%%%%%

DependentDir{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Rum'];
DependentDir{2,1} = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Dis'];
DependentDir{3,1} = [work_dir,'/NodalLevel/Mixed_Effect/HC_Rum'];
DependentDir{4,1} = [work_dir,'/NodalLevel/Mixed_Effect/HC_Dis'];

% only head motion as covariates
OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

% OtherCovariates = [];

% % PALM settings
% PALMSettings.nPerm = 5000;
% PALMSettings.ClusterInference=0;
% PALMSettings.ClusterFormingThreshold=2.3;
% PALMSettings.TFCE=0;
% PALMSettings.FDR=0;
% PALMSettings.TwoTailed=1;
% PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
% PALMSettings.SavePermutations = 1;

mkdir([work_dir, '/NodalLevel/Mixed_Effect/Stats']);
OutputName = [work_dir, '/NodalLevel/Mixed_Effect/Stats/Mixed_Effect_FDR'];

try
    y_MixedEffectsAnalysis_Image(DependentDir,OutputName,'',[],OtherCovariates, []);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end

%%  Paired T MDD
% network level
DependentDir = [];OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/NetworkLevel/Mixed_Effect/MDD_Rum'];
DependentDir{2,1} = [work_dir,'/NetworkLevel/Mixed_Effect/MDD_Dis'];
OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
OutputDir = [work_dir, '/NetworkLevel/PairedT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/PT_MDD'];

try
    y_TTestPaired_Image(DependentDir,OutputName,'',[],OtherCovariates,[]);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end

% Nodal Level
DependentDir = [];OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Rum'];
DependentDir{2,1} = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Dis'];
OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
OutputDir = [work_dir, '/NodalLevel/PairedT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/PT_MDD'];

try
    y_TTestPaired_Image(DependentDir,OutputName,'',[],OtherCovariates, []);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end


% % extract core-MTL functional coonetivity
% FullFC_Rum = y_ReadAll([work_dir,'/NetworkLevel/Mixed_Effect/MDD_Rum']);
% core_MTL_FC_Rum = FullFC_Rum(87,:)';
% FullFC_Dis = y_ReadAll([work_dir,'/NetworkLevel/Mixed_Effect/MDD_Dis']);
% core_MTL_FC_Dis = FullFC_Dis(87,:)';
% core_MTL = [core_MTL_FC_Rum;core_MTL_FC_Dis];
% core_MTL(1:length(core_MTL_FC_Rum),2) = 1;
% core_MTL(length(core_MTL_FC_Rum)+1:length(core_MTL),2) = 2;
% csvwrite([OutputDir,'/core_MTL_FC_MDD.csv'],core_MTL);

%%  Paired T HC
%HC
DependentDir = [];OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/NetworkLevel/Mixed_Effect/HC_Rum'];
DependentDir{2,1} = [work_dir,'/NetworkLevel/Mixed_Effect/HC_Dis'];
OtherCovariates{1,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
OtherCovariates{2,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);
OutputDir =  [work_dir, '/NetworkLevel/PairedT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/PT_HC'];

try
    y_TTestPaired_Image(DependentDir,OutputName,'',[],OtherCovariates, []);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end

% Nodal Level
DependentDir = [];OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/HC_Rum'];
DependentDir{2,1} = [work_dir,'/NodalLevel/Mixed_Effect/HC_Dis'];
OtherCovariates{1,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
OtherCovariates{2,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);
OutputDir = [work_dir, '/NodalLevel/PairedT'];
if ~exist(OutputDir)
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/PT_HC'];

try
    y_TTestPaired_Image(DependentDir,OutputName,'',[],OtherCovariates, []);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end

%% two sample t test rumination
% network level
% PALMSettings.nPerm = 5000;
% PALMSettings.ClusterInference=0;
% PALMSettings.ClusterFormingThreshold=2.3;
% PALMSettings.TFCE=0;
% PALMSettings.FDR=0;
% PALMSettings.TwoTailed=1;
% PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
% PALMSettings.SavePermutations = 1;

DependentDir = [];
OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/NetworkLevel/Mixed_Effect/MDD_Rum'];
DependentDir{2,1} = [work_dir,'/NetworkLevel/Mixed_Effect/HC_Rum'];
OtherCovariates{1,1} = [Age(1:length(find(Dx == 1)),1), Sex(1:length(find(Dx == 1)),1), ...
                        Edu(1:length(find(Dx == 1)),1), HeadMotion(1:length(find(Dx == 1)),1)];
OtherCovariates{2,1} = [Age(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Sex(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Edu(length(find(Dx == 1))+1:length(sub_list),1), ...
                        HeadMotion(length(find(Dx == 1))+1:length(sub_list),1)];

OutputDir = [work_dir, '/NetworkLevel/TwoT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/TwoT_MDD_HC_Rum'];

try
    y_TTest2_Image(DependentDir,OutputName,'',[],OtherCovariates, []);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end

%% Nodal Level
DependentDir = [];
OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Rum'];
DependentDir{2,1} = [work_dir,'/NodalLevel/Mixed_Effect/HC_Rum'];
OtherCovariates{1,1} = [Age(1:length(find(Dx == 1)),1), Sex(1:length(find(Dx == 1)),1), ...
                        Edu(1:length(find(Dx == 1)),1), HeadMotion(1:length(find(Dx == 1)),1)];
OtherCovariates{2,1} = [Age(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Sex(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Edu(length(find(Dx == 1))+1:length(sub_list),1), ...
                        HeadMotion(length(find(Dx == 1))+1:length(sub_list),1)];

OutputDir = [work_dir, '/NodalLevel/TwoT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/TwoT_MDD_HC_Rum'];

try
    y_TTest2_Image(DependentDir,OutputName,'',[],OtherCovariates, []);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end

%% two sample t test distraction
% network level
DependentDir = [];
OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/NetworkLevel/Mixed_Effect/MDD_Dis'];
DependentDir{2,1} = [work_dir,'/NetworkLevel/Mixed_Effect/HC_Dis'];
OtherCovariates{1,1} = [Age(1:length(find(Dx == 1)),1), Sex(1:length(find(Dx == 1)),1), ...
                        Edu(1:length(find(Dx == 1)),1), HeadMotion(1:length(find(Dx == 1)),2)];
OtherCovariates{2,1} = [Age(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Sex(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Edu(length(find(Dx == 1))+1:length(sub_list),1), ...
                        HeadMotion(length(find(Dx == 1))+1:length(sub_list),2)];

OutputDir = [work_dir, '/NetworkLevel/TwoT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/TwoT_MDD_HC_Dis'];

try
    y_TTest2_Image(DependentDir,OutputName,'',[],OtherCovariates, []);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end

% nodal level
DependentDir = [];
OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Dis'];
DependentDir{2,1} = [work_dir,'/NodalLevel/Mixed_Effect/HC_Dis'];
OtherCovariates{1,1} = [Age(1:length(find(Dx == 1)),1), Sex(1:length(find(Dx == 1)),1), ...
                        Edu(1:length(find(Dx == 1)),1), HeadMotion(1:length(find(Dx == 1)),1)];
OtherCovariates{2,1} = [Age(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Sex(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Edu(length(find(Dx == 1))+1:length(sub_list),1), ...
                        HeadMotion(length(find(Dx == 1))+1:length(sub_list),1)];

OutputDir = [work_dir, '/NodalLevel/TwoT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/TwoT_MDD_HC_Dis'];

try
    y_TTest2_Image(DependentDir,OutputName,'',[],OtherCovariates,[]);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end

%% Two sample t test resting state
% network level
DependentDir = [];
OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/NetworkLevel/Mixed_Effect/MDD_Rest'];
DependentDir{2,1} = [work_dir,'/NetworkLevel/Mixed_Effect/HC_Rest'];
OtherCovariates{1,1} = [Age(1:length(find(Dx == 1)),1), Sex(1:length(find(Dx == 1)),1), ...
                        Edu(1:length(find(Dx == 1)),1), HeadMotion(1:length(find(Dx == 1)),2)];
OtherCovariates{2,1} = [Age(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Sex(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Edu(length(find(Dx == 1))+1:length(sub_list),1), ...
                        HeadMotion(length(find(Dx == 1))+1:length(sub_list),2)];

OutputDir = [work_dir, '/NetworkLevel/TwoT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/TwoT_MDD_HC_Rest'];

try
    y_TTest2_Image(DependentDir,OutputName,'',[],OtherCovariates, []);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end

% nodal level
DependentDir = [];
OtherCovariates = [];
DependentDir{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Rest'];
DependentDir{2,1} = [work_dir,'/NodalLevel/Mixed_Effect/HC_Rest'];
OtherCovariates{1,1} = [Age(1:length(find(Dx == 1)),1), Sex(1:length(find(Dx == 1)),1), ...
                        Edu(1:length(find(Dx == 1)),1), HeadMotion(1:length(find(Dx == 1)),1)];
OtherCovariates{2,1} = [Age(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Sex(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Edu(length(find(Dx == 1))+1:length(sub_list),1), ...
                        HeadMotion(length(find(Dx == 1))+1:length(sub_list),1)];

OutputDir = [work_dir, '/NodalLevel/TwoT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/TwoT_MDD_HC_Rest'];

try
    y_TTest2_Image(DependentDir,OutputName,'',[],OtherCovariates,[]);
catch ME
    fprintf(logFileID, '%s - Error occurred: %s\n', datetime, ME.message);
end

%% one sample t test 
% nodal level

% MDD rum
Base = 0;
CovariateDirs = [];
DependentDirs{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Rum'];

OutputDir = [work_dir, '/NodalLevel/OneT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/OneT_MDD_Rum'];

OtherCovariates{1,1} = [Age(1:length(find(Dx == 1)),1), Sex(1:length(find(Dx == 1)),1), ...
                        Edu(1:length(find(Dx == 1)),1), HeadMotion(1:length(find(Dx == 1)),1)];
[TTest1_T,Header] = y_TTest1_Image(DependentDirs,OutputName,[],CovariateDirs,OtherCovariates,Base,[]);

% MDD dis
% nodal level
Base = 0;
CovariateDirs = [];
DependentDirs{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Dis'];

OutputDir = [work_dir, '/NodalLevel/OneT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/OneT_MDD_Dis'];

OtherCovariates{1,1} = [Age(1:length(find(Dx == 1)),1), Sex(1:length(find(Dx == 1)),1), ...
                        Edu(1:length(find(Dx == 1)),1), HeadMotion(1:length(find(Dx == 1)),1)];
[~,~] = y_TTest1_Image(DependentDirs,OutputName,[],CovariateDirs,OtherCovariates,Base,[]);

% HC rum
% nodal level
Base = 0;
CovariateDirs = [];
DependentDirs{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/HC_Rum'];

OutputDir = [work_dir, '/NodalLevel/OneT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/OneT_HC_Rum'];

OtherCovariates{1,1} = [Age(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Sex(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Edu(length(find(Dx == 1))+1:length(sub_list),1), ...
                        HeadMotion(length(find(Dx == 1))+1:length(sub_list),1)];
[~,~] = y_TTest1_Image(DependentDirs,OutputName,[],CovariateDirs,OtherCovariates,Base,[]);

% HC dis
% nodal level
Base = 0;
CovariateDirs = [];
DependentDirs{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/HC_Dis'];

OutputDir = [work_dir, '/NodalLevel/OneT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/OneT_HC_Dis'];

OtherCovariates{1,1} = [Age(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Sex(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Edu(length(find(Dx == 1))+1:length(sub_list),1), ...
                        HeadMotion(length(find(Dx == 1))+1:length(sub_list),1)];
[~,~] = y_TTest1_Image(DependentDirs,OutputName,[],CovariateDirs,OtherCovariates,Base,[]);

% MDD rest
Base = 0;
CovariateDirs = [];
DependentDirs{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/MDD_Rest'];

OutputDir = [work_dir, '/NodalLevel/OneT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/OneT_MDD_Rest'];

OtherCovariates{1,1} = [Age(1:length(find(Dx == 1)),1), Sex(1:length(find(Dx == 1)),1), ...
                        Edu(1:length(find(Dx == 1)),1), HeadMotion(1:length(find(Dx == 1)),1)];
[~,~] = y_TTest1_Image(DependentDirs,OutputName,[],CovariateDirs,OtherCovariates,Base,[]);

% HC rest
Base = 0;
CovariateDirs = [];
DependentDirs{1,1} = [work_dir,'/NodalLevel/Mixed_Effect/HC_Rest'];

OutputDir = [work_dir, '/NodalLevel/OneT'];
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir); 
end
OutputName = [OutputDir,'/OneT_HC_Rest'];

OtherCovariates{1,1} = [Age(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Sex(length(find(Dx == 1))+1:length(sub_list),1), ...
                        Edu(length(find(Dx == 1))+1:length(sub_list),1), ...
                        HeadMotion(length(find(Dx == 1))+1:length(sub_list),1)];
[~,~] = y_TTest1_Image(DependentDirs,OutputName,[],CovariateDirs,OtherCovariates,Base,[]);
