clear all
addpath('./src_classification');
%%
GT_atoms=load('MG_final_atom_coord_inA.mat'); % Ground truth
Type=load('MG_final_atom_types.mat');
Type=Type.final_atom_types;

ref_pos_temp = (GT_atoms.final_atom_coord_inAngstrom); 
ref_pos(1,:)=ref_pos_temp(1,:)-mean(ref_pos_temp(1,:));
ref_pos(2,:)=ref_pos_temp(2,:)-mean(ref_pos_temp(2,:));
ref_pos(3,:)=ref_pos_temp(3,:)-mean(ref_pos_temp(3,:));
ref_pos_all = ref_pos(:,:); % 

ID_Co=find(Type==1);
ID_Pd=find(Type==2);
ID_Pt=find(Type==3);
N=length(ID_Co)
N=length(ID_Pd)
N=length(ID_Pt)

ref_pos_Co = ref_pos(:,ID_Co);% 
ref_pos_Pd = ref_pos(:,ID_Pd);% 
ref_pos_Pt = ref_pos(:,ID_Pt);% 
%%
Model=load('Model_CoPdPt_Dose5p6e4_BM3D_RESIRE.mat'); 

atoms=Model.temp_model;
atoms_type=Model.temp_atomtype;
ID_Co=find(atoms_type==1);
ID_Pd=find(atoms_type==2);% 
ID_Pt=find(atoms_type==3);
%%
atom_pos_all0 = (atoms(:,:)/3-2)*0.333;% 
%%
atom_pos_all(1,:)=atom_pos_all0(2,1:end)-mean(atom_pos_all0(2,1:end));
atom_pos_all(2,:)=atom_pos_all0(1,1:end)-mean(atom_pos_all0(1,1:end));
atom_pos_all(3,:)=-(atom_pos_all0(3,1:end)-mean(atom_pos_all0(3,1:end)));

shift0=[0 0 0]';
atom_pos_all(1,:)=atom_pos_all(1,:)+shift0(1);
atom_pos_all(2,:)=atom_pos_all(2,:)+shift0(2);
atom_pos_all(3,:)=atom_pos_all(3,:)+shift0(3);

atom_pos_Co=atom_pos_all(:,ID_Co); %
atom_pos_Pd=atom_pos_all(:,ID_Pd); % 
atom_pos_Pt=atom_pos_all(:,ID_Pt); % 
%% RMSD of each elements
cutoff=2.0;
[dif_Co,shift_Co,rmsd_Co,Pos_match_Co,~] = Cal_RMSD(ref_pos_Co,atom_pos_Co,cutoff,1);
[dif_Pd,shift_Pd,rmsd_Pd,Pos_match_Pd,~] = Cal_RMSD(ref_pos_Pd,atom_pos_Pd,cutoff,1);
[dif_Pt,shift_Pt,rmsd_Pt,Pos_match_Pt,~] = Cal_RMSD(ref_pos_Pt,atom_pos_Pt,cutoff,1);
%% RMSD of all elements considering atom types
atom_pos_all_classified=[Pos_match_Co Pos_match_Pd Pos_match_Pt];
[dif_all,shift_all,rmsd_all,Pos_match_all,~] = Cal_RMSD(ref_pos_all,atom_pos_all_classified,cutoff,1);
%% RMSD of all elements only consider position as reference
%[~,~] = Cal_RMSD(ref_pos_all,atom_pos_all,cutoff,1);
return
%%



