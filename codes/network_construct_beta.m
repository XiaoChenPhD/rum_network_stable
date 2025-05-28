% Reproduce the GUI work with DPABINet using codes
% First extract ROI signals, then construct network matrix
%
% Xiao Chen
% chenxiaophd@gmail.com
% 250102

%% initialization
clear; clc;

work_dir = 'xxxxxxx';
data_dir = 'xxxxxx';
signal_dir = 'xxxxxxx';
atlas_dir = 'xxxxxxx';

% change this for different conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OutputDir = [signal_dir, '/ROISignals/rest']; 
output_network_dir = [work_dir, '/network_matrix/network_level/rest'];
output_nodal_dir = [work_dir, '/network_matrix/nodal_level/rest'];
data_current_dir = [data_dir, '/Rest/FunSurfWCFS'];
data_currentVolume_dir = [data_dir, '/Rest/FunVoluWCFS'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(OutputDir, 'dir')
    mkdir(OutputDir)
end
if ~exist(output_network_dir, 'dir')
    mkdir(output_network_dir)
end
if ~exist(output_nodal_dir, 'dir')
    mkdir(output_nodal_dir)
end
demog_info_imputed= readtable([signal_dir, '/sample_clean_imputed_v5.csv']);
sub_list = demog_info_imputed.Serial_Number;
StartingDirName = 'FunSurfWCFS';
StartingDirName_Volume = 'FunVoluWCFS';
mask_flag = 'Schaefer400_17networks';

%% extract time series
pool = gcp('nocreate');
if isempty(pool)
    parpool(8);
end
parfor i=1:length(sub_list)
    % left hemisphere
    if ~exist([OutputDir,filesep,'FunSurfLH',filesep, ...
            'ROISignals_',StartingDirName,'_',mask_flag], 'dir')
        mkdir([OutputDir,filesep,'FunSurfLH',filesep, ...
            'ROISignals_',StartingDirName,'_',mask_flag]);
    end
    DirName=dir(fullfile(data_current_dir,['sub-', sub_list{i}], '*fsaverage5_hemi-L*.func.gii'));
    for iFile=1:length(DirName)
        FileName=DirName(iFile).name;
        [ROISignalsSurfLH] = y_ExtractROISignal_Surf( ...
            fullfile(data_current_dir,['sub-', sub_list{i}],FileName), ...
                    {[atlas_dir,filesep,'fsaverage5_lh_Schaefer2018_400Parcels_17Networks_order.label.gii']}, ...
                    [OutputDir,filesep,'FunSurfLH',filesep, ...
                    'ROISignals_',StartingDirName,'_',mask_flag, filesep,sub_list{i}], ...
                    '', 1);
    end
    
    % Right Hemi
    if ~exist([OutputDir,filesep,'FunSurfRH',filesep, 'ROISignals_',StartingDirName,'_',mask_flag], 'dir')
        mkdir([OutputDir,filesep,'FunSurfRH',filesep, ...
            'ROISignals_',StartingDirName,'_',mask_flag]);
    end
    DirName=dir(fullfile(data_current_dir,['sub-', sub_list{i}], '*fsaverage5_hemi-R*.func.gii'));
    for iFile=1:length(DirName)
        FileName=DirName(iFile).name;
        [ROISignalsSurfRH] = y_ExtractROISignal_Surf( ...
            fullfile(data_current_dir,['sub-', sub_list{i}], FileName), ...
                    {[atlas_dir,filesep,'fsaverage5_rh_Schaefer2018_400Parcels_17Networks_order.label.gii']}, ...
                    [OutputDir,filesep,'FunSurfRH',filesep, ...
                    'ROISignals_',StartingDirName,'_',mask_flag, filesep,sub_list{i}], ...
                    '', 1);
    end
    
    % volume
    if ~exist([OutputDir,filesep,'FunVolu',filesep, 'ROISignals_',StartingDirName,'_',mask_flag], 'dir')
       mkdir([OutputDir,filesep,'FunVolu',filesep, 'ROISignals_',StartingDirName,'_',mask_flag]);
    end

   [ROISignalsVolu] = y_ExtractROISignal([data_currentVolume_dir, filesep, 'sub-', sub_list{i}], ...
        {[atlas_dir,filesep,'Tian_Subcortex_S4_3T_2009cAsym.nii']}, ...
        [OutputDir,filesep,'FunVolu',filesep, 'ROISignals_',StartingDirName,'_',mask_flag], ...
        '', 1);
    
    if ~exist([OutputDir, filesep,'ROISignals_SurfLHSurfRHVolu_',StartingDirName], 'dir')
       mkdir([OutputDir, filesep,'ROISignals_SurfLHSurfRHVolu_',StartingDirName]);
    end
    ROISignals = [ROISignalsSurfLH, ROISignalsSurfRH, ROISignalsVolu];
    y_CallSave([OutputDir, filesep,'ROISignals_SurfLHSurfRHVolu_',StartingDirName,filesep, ...
                'ROISignals_',sub_list{i},'.mat'], ROISignals, '');
    y_CallSave([OutputDir, filesep,'ROISignals_SurfLHSurfRHVolu_',StartingDirName,filesep, ...
                'ROISignals_',sub_list{i},'.txt'], ROISignals, ' ''-ASCII'', ''-DOUBLE'',''-TABS''');
end

%% do network construction
pool = gcp('nocreate');
if isempty(pool)
    parpool(4);
end

% nodal level
load([work_dir, '/cfg_TDA_only_surface_nodal.mat']);
DataDir = [OutputDir, filesep,'ROISignals_SurfLHSurfRHVolu_',StartingDirName];
OutDir = output_nodal_dir;
Cfg.SubjectID = sub_list;
[Error, Cfg]=y_NetworkConstruction(Cfg,DataDir,OutDir,[]);

% network level
load([work_dir, '/cfg_TDA_only_surface_network.mat']);
DataDir = [OutputDir, filesep,'ROISignals_SurfLHSurfRHVolu_',StartingDirName];
OutDir = output_network_dir;
Cfg.SubjectID = sub_list;
[Error, Cfg]=y_NetworkConstruction(Cfg,DataDir,OutDir,[]);
