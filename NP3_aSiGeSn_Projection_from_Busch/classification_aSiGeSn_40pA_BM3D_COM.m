close all; clear all;
%%
addpath('./src_classification')
%%
load('Rec_SiGeSn_40pA_10us_BM3D_COM_X4Y4_bgOPT');
%%
load('polyn_tracing_RESIRE_SiGeSn_40pA.mat');
%%
atom_pos = TotPosArr(exitFlagArr==0,:)';
atom_pos_sm = atom_pos./3 - 2;
temp_model = 3*(temp_model_o+2);

[ind1,ind3] = find_commonAtoms(temp_model,temp_atomtype,atom_pos,'ind');
atom_pos(:,ind3) = [];
type1_atom = ones(1,size(atom_pos,2));

tot_atom_pos = cat(2,temp_model,atom_pos);
tot_atom_type = cat(2,temp_atomtype+1,type1_atom);
tot_atom_type(tot_atom_type==4) = 3;

classify_info = struct('Num_species', 3,  'halfSize',  3,  'plothalfSize',  1, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  0,  'separate_part',  100);

[peak_info,intensity_arr,intensity_plot_arr] = ...
    plot_class_hist(FinalVol_single, tot_atom_pos, tot_atom_type, classify_info);

save('peak_info_SiGeSn_40pA_After_BM3D_Alignment.mat',...
    'tot_atom_pos', 'tot_atom_type',...
    'peak_info','intensity_arr','intensity_plot_arr')