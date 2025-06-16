clear all; close all
addpath('./src_classification')
%%
GT_atoms=load('aSiGeSn_8nm_GT.mat'); % Ground truth

Type=(GT_atoms.si_40radius_3Atoms_type)';
ref_pos_temp = (GT_atoms.si_40radius_3Atoms_xyz)'; 

ref_pos(1,:)=ref_pos_temp(2,:)-mean(ref_pos_temp(2,:));
ref_pos(2,:)=ref_pos_temp(1,:)-mean(ref_pos_temp(1,:));
ref_pos(3,:)=ref_pos_temp(3,:)-mean(ref_pos_temp(3,:));

ref_pos_all = ref_pos(:,:); % Cal. RMSD of all atoms.

ID_Si=find(Type==1);
ID_Ge=find(Type==2);
ID_Sn=find(Type==3);
N=length(ID_Si)
N=length(ID_Ge)
N=length(ID_Sn)

ref_pos_Si = ref_pos(:,ID_Si);% 
ref_pos_Ge = ref_pos(:,ID_Ge);% 
ref_pos_Sn = ref_pos(:,ID_Sn);% 
%% load traced results
Model=load('Model_aSiGeSn8nm_FoscCenter_Dose1p6e5_cov35pm_BM3D.mat'); 

atoms=Model.temp_model;
atoms_type=Model.temp_atomtype;
ID_Si=find(atoms_type==1);
ID_Ge=find(atoms_type==2);% 
ID_Sn=find(atoms_type==3);
atom_pos_all = (atoms(:,:)/3-2)*0.347;% 

pos_mean3=[mean(atom_pos_all(1,:)) mean(atom_pos_all(2,:)) mean(atom_pos_all(3,:))];
atom_pos_all(1,:)=atom_pos_all(1,:)-pos_mean3(1);
atom_pos_all(2,:)=atom_pos_all(2,:)-pos_mean3(2);
atom_pos_all(3,:)=atom_pos_all(3,:)-pos_mean3(3);

shift0=[0 0 0]';
atom_pos_all(1,:)=atom_pos_all(1,:)+shift0(1);
atom_pos_all(2,:)=atom_pos_all(2,:)+shift0(2);
atom_pos_all(3,:)=atom_pos_all(3,:)+shift0(3);

atom_pos_Si=atom_pos_all(:,ID_Si); %
atom_pos_Ge=atom_pos_all(:,ID_Ge); % 
atom_pos_Sn=atom_pos_all(:,ID_Sn); % 
%% RMSD of each elements
cutoff=1.5;
[dif_Si,shift_Si,rmsd_Si,Pos_match_Si,~] = Cal_RMSD(ref_pos_Si,atom_pos_Si,cutoff,1);
[dif_Ge,shift_Ge,rmsd_Ge,Pos_match_Ge,~] = Cal_RMSD(ref_pos_Ge,atom_pos_Ge,cutoff,1);
[dif_Sn,shift_Sn,rmsd_Sn,Pos_match_Sn,~] = Cal_RMSD(ref_pos_Sn,atom_pos_Sn,cutoff,1);
%% RMSD of all elements
atom_pos_all_classified=[Pos_match_Si Pos_match_Ge Pos_match_Sn];
[dif_all,shift_all,rmsd_all,Pos_match_all,~] = Cal_RMSD(ref_pos_all,atom_pos_all_classified,cutoff,1);
%%
% save('RMSD_Model_SnGeSn_cov35pmDose1p6e5_S2N2000.mat','dif_Si','dif_Ge','dif_Sn','dif_all','cutoff');


































































return
























return



