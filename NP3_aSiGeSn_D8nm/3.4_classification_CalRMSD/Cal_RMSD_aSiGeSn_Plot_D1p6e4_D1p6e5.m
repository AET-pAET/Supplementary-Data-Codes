clear all; close all;
%%
dif_D1p6e5=load('RMSD_Model_SnGeSn_cov35pmDose1p6e5_S2N2000.mat');
dif_D1p6e4=load('RMSD_Model_SnGeSn_cov35pmDose1p6e4_S2N2000.mat');

dif_AtomTypeAll_D1p6e5=dif_D1p6e5.dif_all;
dif_AtomType1_D1p6e5=dif_D1p6e5.dif_Si;
dif_AtomType2_D1p6e5=dif_D1p6e5.dif_Ge;
dif_AtomType3_D1p6e5=dif_D1p6e5.dif_Sn;

dif_AtomTypeAll_D1p6e4=dif_D1p6e4.dif_all;
dif_AtomType1_D1p6e4=dif_D1p6e4.dif_Si;
dif_AtomType2_D1p6e4=dif_D1p6e4.dif_Ge;
dif_AtomType3_D1p6e4=dif_D1p6e4.dif_Sn;
%%
figure; clf; set(gcf,'Position',[0,0,714,178]);

subplot(1,3,1);
hh=histogram(pdist2(double(dif_AtomType1_D1p6e5'),[0,0,0]), 'BinEdges', 0:0.03:1.0, 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',1.0);
ylabel('Number of atoms')
hold on
hh=histogram(pdist2(double(dif_AtomType1_D1p6e4'),[0,0,0]), 'BinEdges', 0:0.03:1.0, 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',1.0);
xlim([0,1.0])
ylim([0,400])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
lgd = legend('1.6×10^5 e^-/Å^2','1.6×10^4 e^-/Å^2');
legend('boxoff') 
set(legend,'FontSize',9,'linewidth',1.0)
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
lgd.ItemTokenSize = [15, 6] 
set(legend,'Location','NorthWest')

subplot(1,3,2);
hh=histogram(pdist2(double(dif_AtomType2_D1p6e5'),[0,0,0]), 'BinEdges', (0:0.015:0.5), 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',1.0);
xlim([0,0.5])
ylabel('Number of atoms')
hold on
hh=histogram(pdist2(double(dif_AtomType2_D1p6e4'),[0,0,0]), 'BinEdges', (0:0.015:0.5), 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',1.0);
xlim([0,0.5])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
ylim([0,750])

subplot(1,3,3);
hh=histogram(pdist2(double(dif_AtomType3_D1p6e5'),[0,0,0]), 'BinEdges', (0:0.015/2:0.25)/0.25*0.2, 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',1.0);
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
hold on
hh=histogram(pdist2(double(dif_AtomType3_D1p6e4'),[0,0,0]), 'BinEdges', (0:0.015/2:0.25)/0.25*0.2, 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',1.0);
xlim([0,0.20])
ylim([0,550])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);

%%
print('-dbmp', '-r600', 'Fig_2f.bmp');


