% Using codes to get all the graph theory metrics 
%
% Xiao Chen
% 250103
% chenxiaophd@gmail.com

%% initialization
clear; clc;

work_dir = '/mnt/Data3/RfMRILab/ChenX/rumination_network/analyses/redo_DMN';
load('/mnt/Data3/RfMRILab/ChenX/rumination_network/analyses/redo_only_surface/cfg_GTA.mat');
Cfg.ParallelWorkersNumber = 8;

%%
pool = gcp('nocreate');
if isempty(pool)
    parpool(8);
end

DataDir = [work_dir, '/network_matrix/nodal_level/rum'];
OutDir = [work_dir, '/GTA_metrics/rum'];
if ~exist(OutDir, 'dir')
    mkdir(OutDir);
end
[Error, Cfg]=DPABINet_GTA_run(Cfg,DataDir,OutDir,[]);

DataDir = [work_dir, '/network_matrix/nodal_level/dis'];
OutDir = [work_dir, '/GTA_metrics/dis'];
if ~exist(OutDir, 'dir')
    mkdir(OutDir);
end
[Error, Cfg]=DPABINet_GTA_run(Cfg,DataDir,OutDir,[]);