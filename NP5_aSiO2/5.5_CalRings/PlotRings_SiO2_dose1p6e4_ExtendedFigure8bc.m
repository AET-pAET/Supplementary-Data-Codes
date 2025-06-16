clear all; close all
%%
Rings_1p6e4_sl=load('Rings_Model_aSiO2_Dose1p6e4_single_slice_cov35pm_S1p0N1000.mat');
Rings_1p6e4_multi=load('Rings_Model_aSiO2_Dose1p6e4_multislice_cov35pm_S1p0N1000.mat');
%%
cycles_num_ref=Rings_1p6e4_multi.cycles_num_ref;
cycles_sort_ref=Rings_1p6e4_multi.cycles_sort_ref;
Type_all_ref=Rings_1p6e4_multi.Ref_Type;
Pos_all_ref=Rings_1p6e4_multi.ref_pos_all;
%
n=6:1:18;
[GoodRings_No_ref,Atom_type_Ring,Atom_pos_Ring] = Cal_GoodRings(cycles_sort_ref,cycles_num_ref,Pos_all_ref,Type_all_ref,n);
unique(Atom_type_Ring)
%%
cycles_num_1p6e4_multi=Rings_1p6e4_multi.cycles_num_model;
cycles_sort_1p6e4_multi=Rings_1p6e4_multi.cycles_sort_model;
Type_all_1p6e4_multi=Rings_1p6e4_multi.atoms_type_all;
Pos_all_1p6e4_multi=Rings_1p6e4_multi.atom_pos_all;

n=6:1:18;
[GoodRings_No_1p6e4_multi,Atom_type_Ring_1p6e4_multi,Atom_pos_Ring_1p6e4_multi] = Cal_GoodRings(cycles_sort_1p6e4_multi,cycles_num_1p6e4_multi,Pos_all_1p6e4_multi,Type_all_1p6e4_multi,n);
unique(Atom_type_Ring_1p6e4_multi)
%%
cycles_num_1p6e4_sl=Rings_1p6e4_sl.cycles_num_model;
cycles_sort_1p6e4_sl=Rings_1p6e4_sl.cycles_sort_model;
Type_all_1p6e4_sl=Rings_1p6e4_sl.atoms_type_all;
Pos_all_1p6e4_sl=Rings_1p6e4_sl.atom_pos_all;

n=6:1:18;
[GoodRings_No_1p6e4_sl,Atom_type_Ring_1p6e4_sl,Atom_pos_Ring_1p6e4_sl] = Cal_GoodRings(cycles_sort_1p6e4_sl,cycles_num_1p6e4_sl,Pos_all_1p6e4_sl,Type_all_1p6e4_sl,n);
%%
bin=0.25;
width=0.25;
figure; clf; set(gcf,'Position',[0,0,700,150*1]); 
b1 = bar(GoodRings_No_ref(:,1)-bin,GoodRings_No_ref(:,3),'FaceColor',[0 0 160]/255,'EdgeColor',[0 0 0]/255,'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;
hold on

b1 = bar(GoodRings_No_1p6e4_multi(:,1),GoodRings_No_1p6e4_multi(:,3),'FaceColor',[0 102 102]/255,'EdgeColor',[0 0 0]/255,'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;


b1 = bar(GoodRings_No_1p6e4_sl(:,1)+bin,GoodRings_No_1p6e4_sl(:,3),'FaceColor',[237 67 62]/255,'EdgeColor',[0 .0 .0],'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;
ylim([0 1.5])
xlabel(['Number of atoms'])
ylabel('Fraction')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);

legend('GT','Multislice','Single-slice','NumColumns',3);
legend('boxoff') 
set(legend,'FontSize',8,'linewidth',1.0)
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
xticks([6:2:18])
yticks([0:0.3:1.5])
%%
bin=0.17;
width=0.17;
%
[N_ref,X_ref] = hist(cycles_num_ref,6:18);
%
figure; clf; set(gcf,'Position',[0,0,700,150*1]); 
b1 = bar(X_ref-1*bin,N_ref/sum(N_ref),'FaceColor',[0 0 160]/255,'EdgeColor',[0 .0 .0],'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;


[N_multi,X_multi] = hist(cycles_num_1p6e4_multi,6:18);
hold on
b1 = bar(X_multi-0*bin,N_multi/sum(N_multi),'FaceColor',[0 102 102]/255,'EdgeColor',[0 .0 .0],'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;
xlabel(['Number of atoms'])
ylabel('Fraction')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
ylim([0 0.35])

[N_single,X_single] = hist(cycles_num_1p6e4_sl,6:18);
hold on
b1 = bar(X_single+1*bin,N_single/sum(N_single),'FaceColor',[237 67 62]/255,'EdgeColor',[0 .0 .0],'LineWidth',1.0,'FaceAlpha', 0.5);
b1.FaceAlpha = 0.5;
b1.BarWidth = width;
xlabel(['Length of rings'])
ylabel('Fraction')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
ylim([0 0.35])


legend('GT','Multislice','Single-slice','NumColumns',3);

legend('boxoff') 
set(legend,'FontSize',9,'linewidth',1.0)
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
xticks([6:1:18])
xlim([5 19])



