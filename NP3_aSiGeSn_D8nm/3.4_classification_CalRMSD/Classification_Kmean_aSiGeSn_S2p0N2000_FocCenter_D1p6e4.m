clear all; close all
%%
addpath('../3.3_tracing')
addpath('./src_classification');
%% load rec volume
RESIRE=load('Rec_aSiGeSn_FosCetr_Dose1p6e4_cov35pm_BM3D_S2p0N2000_ov3.mat');
FinalVol=RESIRE.FinalVol;
%% load tracing results
load('polyn_tracing_aSiGeSn_FosCetr_Dose1p6e4_cov35pm_BM3D_result.mat') 
atom_pos_raw = TotPosArr(exitFlagArr==0,:); 
%% initial classification of type-2 (Ge) and type-3 (Sn)
classify_info1 = struct('Num_species', 3,  'halfSize',  3,  'plothalfSize',  4, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  120, 'lnorm',2);

[~, temp_atomtype_raw1] = initial_class_kmean_sub(FinalVol, atom_pos_raw', classify_info1);  
%% classification of type-1 (Si) and non-atoms
ID_type2=find(temp_atomtype_raw1==2);
ID_type3=find(temp_atomtype_raw1==3);
ID_type2_type3=[ID_type2 ID_type3];
atom_pos_raw_temp=atom_pos_raw(:,:);
atom_pos_raw_temp(ID_type2_type3,:)=[];
atom_pos_type2_type3=atom_pos_raw(ID_type2_type3,:);

classify_info2 = struct('Num_species', 2,  'halfSize',  1,  'plothalfSize',  4, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  40, 'lnorm',2);
[temp_model2, temp_atomtype2] = initial_class_kmean_sub(FinalVol, atom_pos_raw_temp', classify_info2);

ID_type1=find(temp_atomtype2==1);
ID_type_nonAtoms=ID_type1;
temp_model_nonAtoms=atom_pos_raw_temp(ID_type_nonAtoms,:);
atom_pos_type1=atom_pos_raw_temp(:,:);
atom_pos_type1(ID_type1,:)=[];
%% classification of type-1 (Si), type-2 (Ge) type-3 (Sn)
Atom_pos_all=[atom_pos_type2_type3; atom_pos_type1];

classify_info = struct('Num_species', 3,  'halfSize',  7,  'plothalfSize',  5, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  70, 'lnorm',2);

[temp_model_raw1, temp_atomtype_raw1] = initial_class_kmean_sub(FinalVol, Atom_pos_all', classify_info);

temp_model=[temp_model_raw1];
temp_atomtype=[temp_atomtype_raw1];

temp_model_plot=[temp_model_nonAtoms' temp_model];
temp_atomtype_plot1=temp_atomtype+1;
temp_atomtype_nonAtoms=ID_type_nonAtoms;
temp_atomtype_nonAtoms(:)=1;
temp_atomtype_plot=[temp_atomtype_nonAtoms temp_atomtype_plot1];
%
length(find(temp_atomtype_plot==1))
length(find(temp_atomtype_plot==2))
length(find(temp_atomtype_plot==3))
length(find(temp_atomtype_plot==4))
%% Plot results
classify_info_plot = struct('Num_species', 4,  'halfSize',  3,  'plothalfSize',  4, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  60, 'lnorm',2);

[temp_model_plot, temp_atomtype_plot] = Plot_hist_ColorMap_D1p6e4(FinalVol, temp_model_plot,temp_atomtype_plot, classify_info_plot,0.95,1.1);
%% Output

%save('Model_aSiGeSn8nm_FoscCenter_Dose1p6e4_cov35pm_BM3D.mat',"temp_atomtype","temp_model");






