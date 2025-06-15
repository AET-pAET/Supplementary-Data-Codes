clear all; close all
%%
addpath('./src_classification');
addpath('../5.3_tracing');
%% load rec data
RESIRE=load('Rec_aSiO2_Dose1p6e4_multislice_cov35pm_S1p0N1000_ov3.mat');
FinalVol=RESIRE.FinalVol;
%% load tracing data
load('polyn_tracing_aSiO2_Dose1p6e4_multislice_cov35pm_S1p0N1000_result.mat') %
atom_pos_raw=TotPosArr(exitFlagArr==0,:);
%% initial classification of type-1 (non-Atom), type-2 (O) and type-3 (Si)

classify_info1 = struct('Num_species', 3,  'halfSize',  5,  'plothalfSize',  3, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  120, 'lnorm',2);

[temp_model_classified, temp_atomtype_classified] = initial_class_kmean_sub(FinalVol*100, atom_pos_raw', classify_info1,500);  
%%
temp_model=[temp_model_classified];
temp_atomtype=[temp_atomtype_classified];
return
%% Output
save('Model_aSiO2_Dose1p6e4_multislice_cov35pm_S1p0N1000.mat',"temp_atomtype","temp_model");









