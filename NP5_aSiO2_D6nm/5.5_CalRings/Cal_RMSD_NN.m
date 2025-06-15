function [dif,shift,rmsd_value,Pos_match,count_arr2] = Cal_RMSD_NN(ref_pos,atom_pos,cutoff,plotflag)
Mag_shift=0;
shift=[0 0 0]';

%while Mag_shift>1e-4
m1 = ref_pos - shift;
m2= atom_pos(:,1:end);

difarr = [ ];
count_arr1 = [];
count_arr2 = [];
for i=1:size(m2,2)
    dif=m1-m2(:,i);
    dis=sqrt(sum(dif.^2,1));
    [dis,ind]=min(dis);
    if dis<=cutoff % 0.8 0.45
        difarr=[difarr dif(:,ind)];
        count_arr1 = [count_arr1,ind];
        count_arr2 = [count_arr2,i];
    end
end

shift=shift;%+mean(difarr,2);
Mag_shift=sum(abs(mean(difarr,2)));
%end

% length((count_arr1))/size(m1,2);
% length(unique(count_arr1))/size(m1,2);

fprintf('Number of Atoms = %.0f\n',length((count_arr1)), length(unique(count_arr1)));
fprintf('Accuracy = %.2f\n',length((count_arr1))/size(m1,2)*100, length(unique(count_arr1))/size(m1,2)*100);

dif=difarr;
sig=std(dif,0,2);

rmsd_value = sqrt(sum(sig.^2));
rmsd_value = rmsd_value*100;
fprintf('rmsd = %.2fpm\n',rmsd_value);

Pos_match=atom_pos(:,count_arr2);
Pos_match=Pos_match+shift;
%Pos_match=Pos_match;
if plotflag==1
figure
scatter3(ref_pos(1,:),ref_pos(2,:),ref_pos(3,:),50,'r.')
hold on
%scatter3(atom_pos(1,count_arr2)+shift(1),atom_pos(2,count_arr2)+shift(2),atom_pos(3,count_arr2)+shift(3),20,'k.')
scatter3(Pos_match(1,:),Pos_match(2,:),Pos_match(3,:),20,'k.')

axis equal
axis tight
end
end
    