% Polynomial tracing based on reconstruction of original orientation
clear all; close all; cd(fileparts(mfilename('fullpath')));
addpath('./src');
addpath('../2.3_tracing');
%% load data
RESIRE=load('Rec_aSi9nm_FosCetr_Dose1p6e5_cov35pm_BM3D_S1p0N200_ov3.mat');%
FinalVol=RESIRE.FinalVol;
load('polyn_tracing_aSi9nm_FosCetr_Dose1p6e5_cov35pm_BM3D_result.mat')%reuslts

exitFlagArr0=exitFlagArr(:,1:end); % 
TotPosArr0=TotPosArr(1:end,:);
maxXYZ0=maxXYZ(1:end,:);
atom_pos_tracing = TotPosArr0(exitFlagArr0==0,:); %
%% initial classification
classify_info1 = struct('Num_species', 2,  'halfSize',  5,  'plothalfSize',  3, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  40, 'lnorm',2);

[temp_model_ref, temp_atomtype_ref] = initial_class_kmean_sub_plot_ColorMap(FinalVol*10, atom_pos_tracing', classify_info1,3500);  
%% remove non-atoms
ID_type_nonAtoms=[temp_atomtype_ref==1];
temp_model_Si=atom_pos_tracing(:,:);
temp_model_Si(ID_type_nonAtoms,:)=[];

temp_atomtype_Si=temp_atomtype_ref;
temp_atomtype_Si(ID_type_nonAtoms)=[];

temp_model=[temp_model_Si];
temp_atomtype=[temp_atomtype_Si];

temp_model=temp_model';
temp_atomtype=temp_atomtype';
%% Output
save('Model_aSi9nm_FoscCenter_Dose1p6e5_cov35pm_BM3D.mat',"temp_atomtype","temp_model");
