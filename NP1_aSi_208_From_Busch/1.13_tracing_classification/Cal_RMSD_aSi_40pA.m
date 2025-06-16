clear all; close all;
%%
GT_atoms=load('NP_aSi_208atoms.txt'); % Ground truth
ref_pos_temp = (GT_atoms(:,2:4))';
ref_pos(1,:)=ref_pos_temp(1,:)-mean(ref_pos_temp(1,:));
ref_pos(2,:)=ref_pos_temp(2,:)-mean(ref_pos_temp(2,:));
ref_pos(3,:)=ref_pos_temp(3,:)-mean(ref_pos_temp(3,:));
%% Classified atoms
Model=load('Model_aSi_40pA_FromBM3DRESIRE_check.mat');  
atoms=Model.temp_model;
atoms_type=Model.temp_atomtype;

atom_pos_all = (atoms/3-2)*0.195;
atom_pos_temp=atom_pos_all(:,1:end);
atom_pos(1,:)=atom_pos_temp(2,:)-mean(atom_pos_temp(2,:));
atom_pos(2,:)=(atom_pos_temp(1,:)-mean(atom_pos_temp(1,:)));
atom_pos(3,:)=-(atom_pos_temp(3,:)-mean(atom_pos_temp(3,:)));

shift0=[0 0 0]';
atom_pos(1,:)=atom_pos(1,:)+shift0(1);
atom_pos(2,:)=atom_pos(2,:)+shift0(2);
atom_pos(3,:)=atom_pos(3,:)+shift0(3);
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
dif=difarr;
sig=std(dif,0,2);
rmsd_value = sqrt(sum(sig.^2));
fprintf('rmsd = %.2fpm\n',rmsd_value*100);
%%
figure
scatter3(ref_pos(1,:),ref_pos(2,:),ref_pos(3,:),100,'b.')
hold on
scatter3(atom_pos(1,count_arr2)+shift(1),atom_pos(2,count_arr2)+shift(2),atom_pos(3,count_arr2)+shift(3),200,'r.')
axis equal
axis tight
xlim([-10 10])
ylim([-10 10])
zlim([-10 10])
%% Plot ExtDataFig.3e
figure(11); clf; set(gcf,'Position',[0,100,350,220]);
hh = histogram(pdist2(double(dif'),[0,0,0]),0:0.01:0.15,'linewidth',1);
hh.FaceColor = [0 0 160]/255;
hh.FaceAlpha = 0.5; 
xlim([0,0.15])
ylim([0 50])
xlabel(['Deviation(',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
%%
% for i=1:208
% d(i)=sqrt(dif(1,i)^2+dif(2,i)^2+dif(3,i)^2);
% end
% 
% mean(d)