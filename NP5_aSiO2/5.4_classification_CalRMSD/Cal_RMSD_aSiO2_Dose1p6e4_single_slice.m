clear all
%%
GT_atoms=load('model_GT.mat'); % Ground truth
Ref_Type=GT_atoms.atomic_number;

ref_pos = (GT_atoms.pos)'; 
ref_pos0(1,:)=ref_pos(2,:)-mean(ref_pos(2,:));
ref_pos0(2,:)=ref_pos(1,:)-mean(ref_pos(1,:));
ref_pos0(3,:)=ref_pos(3,:)-mean(ref_pos(3,:));

ref_pos_all = ref_pos0(:,:); 

ID_O=find(Ref_Type==8);
ID_Si=find(Ref_Type==14);
N=length(ID_O)
N=length(ID_Si)

ref_pos_O = ref_pos0(:,ID_O);% 
ref_pos_Si = ref_pos0(:,ID_Si);% 
%% Classified atoms
Model=load('Model_aSiO2_Dose1p6e4_single_slice_cov35pm_S1p0N1000.mat'); 

atoms=Model.temp_model;
atoms_type=Model.temp_atomtype;
ID_nonAtom=find(atoms_type==1);
ID_O=find(atoms_type==2);% 
ID_Si=find(atoms_type==3);
ID_all=[ID_O ID_Si];
%%
atom_pos_all = (atoms(:,ID_all)/3-2)*0.34517;% 
atoms_type_all=atoms_type(ID_all);
pos_mean3=[mean(atom_pos_all(1,:)) mean(atom_pos_all(2,:)) mean(atom_pos_all(3,:))];
atom_pos_all(1,:)=atom_pos_all(1,:)-pos_mean3(1);
atom_pos_all(2,:)=atom_pos_all(2,:)-pos_mean3(2);
atom_pos_all(3,:)=atom_pos_all(3,:)-pos_mean3(3);

shift0=[0 0 0]';
atom_pos_all(1,:)=atom_pos_all(1,:)+shift0(1);
atom_pos_all(2,:)=atom_pos_all(2,:)+shift0(2);
atom_pos_all(3,:)=atom_pos_all(3,:)+shift0(3);

atom_pos_O=atom_pos_all(:,1:length(ID_O)); %
atom_pos_Si=atom_pos_all(:,length(ID_O)+1:end); % 
%% RMSD of each elements
cutoff=1.0;
[dif_Si,shift_Si,rmsd_Si,Pos_match_Si,~] = Cal_RMSD(ref_pos_Si,atom_pos_Si,cutoff,1);
[dif_O,shift_O,rmsd_O,Pos_match_O,~] = Cal_RMSD(ref_pos_O,atom_pos_O,cutoff,1);
%% RMSD of all elements considering atom types
atom_pos_all_classified=[Pos_match_Si Pos_match_O];
[dif_all,shift_all,rmsd_all,Pos_match_all,~] = Cal_RMSD(ref_pos_all,atom_pos_all_classified,cutoff,1);
%% RMSD of all elements only consider position as reference
[~,~] = Cal_RMSD(ref_pos_all,atom_pos_all,cutoff,1);
%%
return




