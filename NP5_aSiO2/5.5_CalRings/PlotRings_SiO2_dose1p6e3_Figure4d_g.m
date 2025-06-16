clear all; close all;
%%
Rings_1p6e3_sl=load('Rings_Model_aSiO2_Dose1p6e3_single_slice_cov35pm_S2p0N1000.mat');
Rings_1p6e3_multi=load('Rings_Model_aSiO2_Dose1p6e3_multislice_cov35pm_S1p0N1000.mat');
%%
cycles_num_ref=Rings_1p6e3_multi.cycles_num_ref;
cycles_sort_ref=Rings_1p6e3_multi.cycles_sort_ref;
Type_all_ref=Rings_1p6e3_multi.Ref_Type;
Pos_all_ref=Rings_1p6e3_multi.ref_pos_all;
%
n=6:1:18;
[GoodRings_No_ref,Atom_type_Ring,Atom_pos_Ring] = Cal_GoodRings(cycles_sort_ref,cycles_num_ref,Pos_all_ref,Type_all_ref,n);
unique(Atom_type_Ring)
%% Plot Rings
% n0=1-17, id=12; n0=18, id=13
n0=6;
id=12;
[GoodRings_No_ref0,Atom_type_Ring0,Atom_pos_Ring0] = Cal_GoodRings_Plot(cycles_sort_ref,cycles_num_ref,Pos_all_ref,Type_all_ref,n0,id); % ~12
unique(Atom_type_Ring0)
Cal_Rings_Plot(Atom_type_Ring0',Atom_pos_Ring0',[8 14],6,90,0);
%%
cycles_num_1p6e3_multi=Rings_1p6e3_multi.cycles_num_model;
cycles_sort_1p6e3_multi=Rings_1p6e3_multi.cycles_sort_model;
Type_all_1p6e3_multi=Rings_1p6e3_multi.atoms_type_all;
Pos_all_1p6e3_multi=Rings_1p6e3_multi.atom_pos_all;

n=6:1:18;
[GoodRings_No_1p6e3_multi,Atom_type_Ring_1p6e3_multi,Atom_pos_Ring_1p6e3_multi] = Cal_GoodRings(cycles_sort_1p6e3_multi,cycles_num_1p6e3_multi,Pos_all_1p6e3_multi,Type_all_1p6e3_multi,n);
unique(Atom_type_Ring_1p6e3_multi)
%%
cutoff=0.5;
plotflag=0;
for i=1:n0
[dif,shift,rmsd_value,Ring_match_1p6e3_multi_temp,ID_arr_1p6e3_multi_temp] = Cal_RMSD_NN(Atom_pos_Ring0(:,i),Pos_all_1p6e3_multi,cutoff,plotflag);
%
Ring_match_1p6e3_multi(:,i)=Ring_match_1p6e3_multi_temp;
Type_match_1p6e3_multi(i)=Type_all_1p6e3_multi(ID_arr_1p6e3_multi_temp);
end

unique(Type_match_1p6e3_multi)
%%
close all;
figure('color','white');
set(gcf,'Position',[100,100,210*1,210*2]);
subplot(2,1,1);
Cal_Rings_Plot(Atom_type_Ring0',Atom_pos_Ring0',[8 14],6.5,11,35);
subplot(2,1,2);
Cal_Rings_Plot(Type_match_1p6e3_multi',Ring_match_1p6e3_multi',[2 3],6.5,11,35);
%print('-djpeg', '-r1200', 'Fig_4e_marker.jpg');
%%
% N=6, view(11,35)
% N=8, view(11,70)
% N=10, view(32,150)
% N=12, view(134,77)
% N=14, view(135,50)
% N=16, view(204,20)
% N=18, view(-22,27)
%%
cycles_num_1p6e3_sl=Rings_1p6e3_sl.cycles_num_model;
cycles_sort_1p6e3_sl=Rings_1p6e3_sl.cycles_sort_model;
Type_all_1p6e3_sl=Rings_1p6e3_sl.atoms_type_all;
Pos_all_1p6e3_sl=Rings_1p6e3_sl.atom_pos_all;

n=6:1:18;
[GoodRings_No_1p6e3_sl,Atom_type_Ring_1p6e3_sl,Atom_pos_Ring_1p6e3_sl] = Cal_GoodRings(cycles_sort_1p6e3_sl,cycles_num_1p6e3_sl,Pos_all_1p6e3_sl,Type_all_1p6e3_sl,n);
%%
bin=0.25;
width=0.25;
figure(11); clf; set(gcf,'Position',[0,0,380,150*1]); % 450
b1 = bar(GoodRings_No_ref(:,1)-bin,GoodRings_No_ref(:,3),'FaceColor',[0 0 160]/255,'EdgeColor',[0 0 0]/255,'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;
hold on

b1 = bar(GoodRings_No_1p6e3_multi(:,1),GoodRings_No_1p6e3_multi(:,3),'FaceColor',[0 102 102]/255,'EdgeColor',[0 0 0]/255,'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;

b1 = bar(GoodRings_No_1p6e3_sl(:,1)+bin,GoodRings_No_1p6e3_sl(:,3),'FaceColor',[237 67 62]/255,'EdgeColor',[0 .0 .0],'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;
ylim([0 1.5])
xlabel(['Number of atoms'])
ylabel('Fraction')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);

lgd=legend('GT','Multislice','Single-slice','NumColumns',3);%  without probe

legend('boxoff') 
lgd.ItemTokenSize = [20, 6]

set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
%set(legend,'Location','NorthWest') 
xticks([6:2:18])
yticks([0:0.5:1.5])
% print('-dbmp', '-r1200', 'Fig_4g.bmp');
% print('-djpeg', '-r2400', 'Fig_4g.jpg');
%%
bin=0.25;
width=0.25;
%
[N_ref,X_ref] = hist(cycles_num_ref,6:18);
%
figure(11); clf; set(gcf,'Position',[0,0,380,150*1]); 
b1 = bar(X_ref-bin,N_ref/sum(N_ref),'FaceColor',[0 0 160]/255,'EdgeColor',[0 .0 .0],'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;

%
[N_multi,X_multi] = hist(cycles_num_1p6e3_multi,6:18);
hold on
b1 = bar(X_multi,N_multi/sum(N_multi),'FaceColor',[0 102 102]/255,'EdgeColor',[0 .0 .0],'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;
xlabel(['Number of atoms'])
ylabel('Fraction')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
ylim([0 0.38])
%
[N_single,X_single] = hist(cycles_num_1p6e3_sl,6:18);
hold on
b1 = bar(X_single+bin,N_single/sum(N_single),'FaceColor',[237 67 62]/255,'EdgeColor',[0 .0 .0],'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;
xlabel(['Number of atoms'])
ylabel('Fraction')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
ylim([0 0.38])

lgd=legend('GT','Multislice','Single-slice','NumColumns',3);%  without probe


legend('boxoff') 
lgd.ItemTokenSize = [20, 6]
set(legend,'FontSize',9,'linewidth',1.0)
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(legend,'Location','NorthWest') 
xticks([6:1:18])
xlim([5 19])
%
% print('-dbmp', '-r1200', 'Fig_4f.bmp');
% print('-djpeg', '-r2400', 'Fig_4f.jpg');
