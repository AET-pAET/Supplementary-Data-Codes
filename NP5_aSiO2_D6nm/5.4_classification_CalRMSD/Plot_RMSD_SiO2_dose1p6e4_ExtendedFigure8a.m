clear all; close all;
%%
RMSD_D1p6e4_cov35pm_multi=load('RMSD_Model_aSiO2_Dose1p6e4_multislice_cov35pm_S1p0N1000.mat');
RMSD_D1p6e4_cov35pm_multi_Si=RMSD_D1p6e4_cov35pm_multi.dif_Si;
RMSD_D1p6e4_cov35pm_multi_O=RMSD_D1p6e4_cov35pm_multi.dif_O;
RMSD_D1p6e4_cov35pm_multi_All=RMSD_D1p6e4_cov35pm_multi.dif_all;
%%
RMSD_D1p6e4_cov35pm_single=load('RMSD_Model_aSiO2_Dose1p6e4_single_slice_cov35pm_S1p0N1000.mat');
RMSD_D1p6e4_cov35pm_single_Si=RMSD_D1p6e4_cov35pm_single.dif_Si;
RMSD_D1p6e4_cov35pm_single_O=RMSD_D1p6e4_cov35pm_single.dif_O;
RMSD_D1p6e4_cov35pm_single_All=RMSD_D1p6e4_cov35pm_single.dif_all;
%%

figure; clf; set(gcf,'Position',[0,0,230*3.5*0.9,210*0.9]); 

subplot(1,3,2);
hh=histogram(pdist2(double(RMSD_D1p6e4_cov35pm_multi_O'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.5*0.8, 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
ylabel('Number of atoms')
hold on
hh=histogram(pdist2(double(RMSD_D1p6e4_cov35pm_single_O'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.5*0.8, 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
xlim([0,0.4])
ylim([0,900])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);

subplot(1,3,1);
hh=histogram(pdist2(double(RMSD_D1p6e4_cov35pm_multi_Si'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.5*0.8/0.4*0.3, 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
ylabel('Number of atoms')
hold on
hh=histogram(pdist2(double(RMSD_D1p6e4_cov35pm_single_Si'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.5*0.8/0.4*0.3, 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
xlim([0,0.3])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
ylim([0,600])

leg=legend('Multislice','Single-slice');
legend('boxoff') 
set(legend,'FontSize',8,'linewidth',1.0)
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(legend,'Location','NorthWest') 
leg.ItemTokenSize = [18, 6];

subplot(1,3,3);
hh=histogram(pdist2(double(RMSD_D1p6e4_cov35pm_multi_All'),[0,0,0]), 'BinEdges', (0:0.015:0.5)*0.8, 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
ylabel('Number of atoms')
hold on
hh=histogram(pdist2(double(RMSD_D1p6e4_cov35pm_single_All'),[0,0,0]), 'BinEdges', (0:0.015:0.5)*0.8, 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
xlim([0,0.4])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);



