clear all
cd(fileparts(mfilename('fullpath')));
addpath('./');
%%
GT_atoms=load('aSi9nm_model_GT.txt'); % Ground truth
N=length(GT_atoms(:,1))
ref_pos1 = (GT_atoms(:,2:4))';
ref_pos(1,:)=ref_pos1(2,:)-mean(ref_pos1(2,:));
ref_pos(2,:)=ref_pos1(1,:)-mean(ref_pos1(1,:));
ref_pos(3,:)=ref_pos1(3,:)-mean(ref_pos1(3,:));
%% Classified atoms
Model=load('Model_aSi9nm_FoscCenter_Dose1p6e5_cov35pm_BM3D.mat'); 

atoms=Model.temp_model;
atoms_type=Model.temp_atomtype;
ID_temp=find(atoms_type==2);

atom_pos_all = (atoms/3-2)*0.347;% 
atom_pos=atom_pos_all(:,ID_temp); 
pos_mean3=[mean(atom_pos(1,:)) mean(atom_pos(2,:)) mean(atom_pos(3,:))];

atom_pos(1,:)=atom_pos(1,:)-pos_mean3(1);
atom_pos(2,:)=atom_pos(2,:)-pos_mean3(2);
atom_pos(3,:)=atom_pos(3,:)-pos_mean3(3);

shift0=[0 0 0]';
atom_pos(1,:)=atom_pos(1,:)+shift0(1);
atom_pos(2,:)=atom_pos(2,:)+shift0(2);
atom_pos(3,:)=atom_pos(3,:)+shift0(3);
%return
%%
Mag_shift=1;
shift=[0 0 0]';

while Mag_shift>1e-4
m1 = ref_pos - shift;
m2= atom_pos(:,1:end);

difarr = [ ];
count_arr1 = [];
count_arr2 = [];
for i=1:size(m2,2)
    dif=m1-m2(:,i);
    dis=sqrt(sum(dif.^2,1));
    [dis,ind]=min(dis);
    if dis<=1.0
        difarr=[difarr dif(:,ind)];
        count_arr1 = [count_arr1,ind];
        count_arr2 = [count_arr2,i];
    end
end

shift=shift+mean(difarr,2);
Mag_shift=sum(abs(mean(difarr,2)));
end
disp(num2str(numel(count_arr1)/size(m1,2)))
%%
dif=difarr;
sig=std(dif,0,2);

rmsd_value = sqrt(sum(sig.^2));
fprintf('rmsd = %.2fpm\n',rmsd_value*100);
%% Output 
% save('dif_Dose1p6e5_AtomTypeAll.mat','dif')






