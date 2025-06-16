clear all; close all
%%
addpath('./src');
addpath('./src_classification');
%% load rec data
RESIRE=load('Rec_CoPdPt_NoRadmDefocus_BM3D_ExpAngle_S1p0N1000_ov3.mat');
FinalVol=RESIRE.FinalVol;
%%
load('polyn_tracing_CoPdPt_NoRadmDefocus_BM3D_ExpAngle_S1p0N1000_result.mat') % 
atom_pos_raw=TotPosArr(exitFlagArr==0,:);% 
%%  initial classification
classify_info1 = struct('Num_species', 3,  'halfSize',  3,  'plothalfSize',  4, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  120, 'lnorm',2);

[~, temp_atomtype_raw] = initial_class_kmean_sub(FinalVol, atom_pos_raw', classify_info1);  
%% 
ID_type2=find(temp_atomtype_raw==2);
ID_type3=find(temp_atomtype_raw==3);
ID_type2_type3=[ID_type2 ID_type3];

atom_pos_type1_nonAtoms=atom_pos_raw(:,:);
atom_pos_type1_nonAtoms(ID_type2_type3,:)=[];
atom_pos_type2_type3=atom_pos_raw(ID_type2_type3,:);

classify_info2 = struct('Num_species', 2,  'halfSize',  3,  'plothalfSize',  4, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  40, 'lnorm',2);
 
[~, temp_atomtype_type1_nonAtoms] = initial_class_kmean_sub(FinalVol, atom_pos_type1_nonAtoms', classify_info2);

ID_type_nonAtoms=find(temp_atomtype_type1_nonAtoms==1);
temp_model_nonAtoms=atom_pos_type1_nonAtoms(ID_type_nonAtoms,:);
atom_pos_type1=atom_pos_type1_nonAtoms(:,:);
atom_pos_type1(ID_type_nonAtoms,:)=[];
%% classification of type-1 (Co), type-2 (Pd) type-3 (Pt)
Atom_pos_all=[atom_pos_type2_type3; atom_pos_type1];

classify_info = struct('Num_species', 3,  'halfSize',  9,  'plothalfSize',  3, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  70, 'lnorm',2);

[temp_model_classified, temp_atomtype_classified] = initial_class_kmean_sub(FinalVol, Atom_pos_all', classify_info);
%% Output
temp_model=[temp_model_classified];
temp_atomtype=[temp_atomtype_classified];
%% plot hist
temp_model_plot=[temp_model_nonAtoms' temp_model];
temp_atomtype_plot1=temp_atomtype+1;
temp_atomtype_nonAtoms=ID_type_nonAtoms;
temp_atomtype_nonAtoms(:)=1;
temp_atomtype_plot=[temp_atomtype_nonAtoms temp_atomtype_plot1];

length(find(temp_atomtype_plot==1))
length(find(temp_atomtype_plot==2))
length(find(temp_atomtype_plot==3))
length(find(temp_atomtype_plot==4))
%
classify_info_plot = struct('Num_species', 4,  'halfSize',  3,  'plothalfSize',  4, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  100, 'lnorm',2);

[temp_model_plot, temp_atomtype_plot] = Plot_hist_CoPdPt_ExtFig7c(FinalVol, temp_model_plot,temp_atomtype_plot, classify_info_plot,1500,18);
return
%% Output
save('Model_CoPdPt_NoRadmDefocus_BM3D_ExpAngle_S1p0N1000.mat',"temp_atomtype","temp_model");






