clear all; close all
%%
addpath('./src');
addpath('../1.3_tracing');
%% load rec volume 
RESIRE=load('Rec_aSi2nm_FosCetr_Dose1p6e4_cov35pm_BM3D_S0p25N200_ov3.mat');
FinalVol=RESIRE.FinalVol;
%% load tracing data
load('polyn_tracing_aSi2nm_FosCetr_Dose1p6e4_cov35pm_BM3D.mat') %

exitFlagArr0=exitFlagArr(:,1:end); % 
TotPosArr0=TotPosArr(1:end,:);
maxXYZ0=maxXYZ(1:end,:);
atom_pos_tracing = TotPosArr0(exitFlagArr0==0,:); %
%% classification
classify_info1 = struct('Num_species', 2,  'halfSize',  3,  'plothalfSize',  3, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  40, 'lnorm',2);

[temp_model_ref, temp_atomtype_ref] = initial_class_kmean_sub_plot_ColorMap(FinalVol*10, atom_pos_tracing', classify_info1);
%% speration of type-1 (Si) and non-atoms
ID_type_nonAtoms=[temp_atomtype_ref==1];
atom_pos_Si=atom_pos_tracing(:,:);
atom_pos_Si(ID_type_nonAtoms,:)=[];

classify_info2 = struct('Num_species', 1,  'halfSize',  5,  'plothalfSize',  5, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  20, 'lnorm',2);

[temp_model_Si, temp_atomtype_Si] = initial_class_kmean_sub(FinalVol*10, atom_pos_Si', classify_info2);

temp_model=[temp_model_Si];
temp_atomtype=[temp_atomtype_Si];
%% Output
save('Model_aSi2nm_FoscCenter_Dose1p6e4_cov35pm_BM3D.mat',"temp_atomtype","temp_model");








