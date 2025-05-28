function [corrected_p_matrix, adjusted_p_matrix] = c_fdr_network(NetworkMatrix, FDRQ, StatOpt, flag_corrected)
% function FinalCorrectedMatrix = c_fdr_network(NetworkMatrix, FDRQ, StatOpt)
% do fdr correction according to the setted FDRQ on the output network
% stats matrix from dpabiNet
% Inputs:
% NetworkMatrix: the stats matrix from dpabiNetwork, containing t values
% FDRQ: the threshold of the fdr correction
% StatOpt: a structure from the dpabiNet, denoting the statistic
%          properties. It must consist of a field named "TailedFlag". 
%          1, one-tailed, 2, two-tailed
% flag_corrected: 1: network level; 2: node level.
%
% Outputs:
% corrected_p_matrix: a matrix only consisting of 1s and 0s, 1s denoting
%                       edges that survive the fdr correction
% 
% adjusted_p_matrix: the FDR adjusted p values
%
% Xiao Chen, PhD
% 250120
%
% V2: use mafdr to do FDR correction to get the adjusted p values and fix a
% potential bug
% 
% Xiao Chen
% 250121
% chenxiaophd@gmail.com

if flag_corrected == 1
    Vec = tril(NetworkMatrix);
    VecP = w_StatToP(Vec, StatOpt); %Convert the Stat values to p values
    VecP = VecP(tril(ones(size(VecP))) == 1); % to a vector
    adjusted_VecP = mafdr(VecP, 'BHFDR', true);
    % binarized corrected map
    corrected_VecP = adjusted_VecP;
    corrected_VecP(corrected_VecP > FDRQ) = 0;
    corrected_VecP(corrected_VecP~=0) = 1;
    % map back to the matrix
    lower_tri_indices = tril(true(size(NetworkMatrix)));
    adjusted_p_matrix = zeros(size(NetworkMatrix));
    adjusted_p_matrix(lower_tri_indices) = adjusted_VecP;
    adjusted_p_matrix = adjusted_p_matrix + adjusted_p_matrix' - diag(diag(adjusted_p_matrix));
    % binarized
    corrected_p_matrix = zeros(size(NetworkMatrix));
    corrected_p_matrix(lower_tri_indices) = corrected_VecP;
    corrected_p_matrix = corrected_p_matrix + corrected_p_matrix' - diag(diag(corrected_p_matrix));
elseif flag_corrected == 2
    Vec = tril(NetworkMatrix, -1);
    VecP = w_StatToP(Vec, StatOpt); %Convert the Stat values to p values
    VecP = VecP(tril(ones(size(VecP)), -1) == 1); % to a vector
    adjusted_VecP = mafdr(VecP, 'BHFDR', true);
    % binarized corrected map
    corrected_VecP = adjusted_VecP;
    corrected_VecP(corrected_VecP > FDRQ) = 0;
    corrected_VecP(corrected_VecP~=0) = 1;
    % map back to the matrix
    lower_tri_indices = tril(true(size(NetworkMatrix)), -1);
    adjusted_p_matrix = zeros(size(NetworkMatrix));
    adjusted_p_matrix(lower_tri_indices) = adjusted_VecP;
    adjusted_p_matrix = adjusted_p_matrix + adjusted_p_matrix' - diag(diag(adjusted_p_matrix));
    % setting all the diagnal elements as 1s
    adjusted_p_matrix(1:size(adjusted_p_matrix,1)+1:end) = 1;
    % binarized
    corrected_p_matrix = zeros(size(NetworkMatrix));
    corrected_p_matrix(lower_tri_indices) = corrected_VecP;
    corrected_p_matrix = corrected_p_matrix + corrected_p_matrix' - diag(diag(corrected_p_matrix));
end