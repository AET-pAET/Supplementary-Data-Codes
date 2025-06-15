% Polynomial tracing based on reconstruction of original orientation
clear all; close all; cd(fileparts(mfilename('fullpath')));
%%
addpath('./src');
%% load rec data
RESIRE=load('Rec_Si_40pA_BM3D_ov3.mat');
FinalVol=RESIRE.FinalVol;
%% load tracing data
load('polyn_tracing_40pA_N5000_result.mat') %
atom_pos_raw=TotPosArr(exitFlagArr==0,:);
%% classification of NonAtoms and Si
classify_info1 = struct('Num_species', 2,  'halfSize',  3,  'plothalfSize',  3, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  20, 'lnorm',2);
 
[~, temp_atomtype_raw] = plot_ExtDataFig3d(FinalVol*10, atom_pos_raw', classify_info1);  
%% remove non-atoms
ID_type_NonAtoms=find(temp_atomtype_raw==1);
atom_pos_Si=atom_pos_raw(:,:);
atom_pos_Si(ID_type_NonAtoms,:)=[];

classify_info2 = struct('Num_species', 1,  'halfSize',  5,  'plothalfSize',  5, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  30, 'lnorm',2);

[temp_model_classified, temp_atomtype_classified] = initial_class_kmean_sub(FinalVol, atom_pos_Si', classify_info2);
%% Output
temp_model=[temp_model_classified];
temp_atomtype=[temp_atomtype_classified];
save('Model_aSi_40pA_FromBM3DRESIRE_check.mat',"temp_atomtype","temp_model");










