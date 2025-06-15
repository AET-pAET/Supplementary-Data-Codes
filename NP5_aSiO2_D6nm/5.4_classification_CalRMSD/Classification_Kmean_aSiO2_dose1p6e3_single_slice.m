% Polynomial tracing based on reconstruction of original orientation
clear all; close all
%%
addpath('./src_classification');
addpath('../5.3_tracing');
%% load data
RESIRE=load('Rec_aSiO2_Dose1p6e3_single_slice_cov35pm_S2p0N1000_ov3.mat'); 
FinalVol=RESIRE.FinalVol;
%%
load('polyn_tracing_aSiO2_Dose1p6e3_single_slice_cov35pm_S2p0N1000_result.mat');
atom_pos_raw = TotPosArr(exitFlagArr==0,:); %
%% initial classification of non-Atoms (type-1) and atoms.

classify_info1 = struct('Num_species', 2,  'halfSize',  3,  'plothalfSize',  2, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  120, 'lnorm',2);

[temp_model_raw, temp_atomtype_raw] = initial_class_kmean_sub(FinalVol*100, atom_pos_raw', classify_info1,250);  
%% remove nonAtoms and speration of Si and O
ID_type_nonAtoms=find(temp_atomtype_raw==1);
atom_pos_temp=atom_pos_raw(:,:);
atom_pos_temp(ID_type_nonAtoms,:)=[];

classify_info2 = struct('Num_species', 2,  'halfSize',  3,  'plothalfSize',  2, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  70, 'lnorm',2);
% 
 [temp_model_classified, temp_atomtype_classified] = initial_class_kmean_sub(FinalVol*100, atom_pos_temp', classify_info2,250);
%% Plot histgram
ID_nonAtom=find(temp_atomtype_raw==1);
temp_atomtype_nonAtom=temp_atomtype_raw(ID_nonAtom);
Temp_model_nonAtom=temp_model_raw(:,ID_nonAtom);
temp_atomtype2_plot=temp_atomtype_classified+1;
temp_atomtype_plot=[temp_atomtype_nonAtom temp_atomtype2_plot];
Temp_model_plot=[Temp_model_nonAtom temp_model_classified];

classify_info_plot = struct('Num_species', 3,  'halfSize',  3,  'plothalfSize',  2, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  100, 'lnorm',2);

[~, ~] = Hist_plot_aSiO2_Dose1p6e3_single_slice(FinalVol*1000, Temp_model_plot, temp_atomtype_plot, classify_info_plot,600,12);
%%
temp_model=[temp_model_classified];
temp_atomtype=[temp_atomtype_classified];
return
%% Output
save('Model_aSiO2_Dose1p6e3_single_slice_cov35pm_S2p0N1000.mat',"temp_atomtype","temp_model");







